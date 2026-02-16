from __future__ import annotations

import glob
import os
import subprocess
from typing import Any, Dict, List, Optional, Tuple, Union
from pathlib import Path

import yaml

from frogger.common.codon_table import load_codon_usage_xlsx
from frogger.common.gg.splitters import split_cds_to_fragment_dicts
from frogger.common.io import read_fasta_records, write_fasta_records, write_tsv
from frogger.gg_adapters import build_cloning_sequences_for_fragments
from frogger.overhangs import load_pairing_matrix_xlsx

# Converge on SHUFFLE GG logic (overhang selection + primers)
from frogger.common.gg.packager import (
    build_primers_by_pos_for_fragments,
    compute_overhangs_by_frag_pos,
    select_global_overhangs_and_forced_codons,
)

from .mutate import (
    Variant,
    generate_explicit_variants,
    generate_saturation_variants,
    translate_cds,
)

MutationSpec = Tuple[str, Union[List[int], List[str]]]


def load_config(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise ValueError("Config must be a YAML mapping.")
    return cfg


def _merge_cut_points(cfg: Dict[str, Any], gene_id: str) -> Tuple[List[int], bool]:
    splits = cfg.get("splits", {}) or {}
    global_cps = list(splits.get("global_cut_points", []) or [])
    per_gene = (splits.get("per_gene_cut_points", {}) or {})
    enforce = bool(splits.get("enforce_frame", True))
    if gene_id in per_gene and per_gene[gene_id] is not None:
        return list(per_gene[gene_id]), enforce
    return global_cps, enforce


def _collect_forbidden_motifs(cfg: Dict[str, Any]) -> List[str]:
    codonopt_cfg = cfg.get("codonopt", {}) or {}
    avoid = codonopt_cfg.get("avoid_motifs", []) or []
    final_forbidden = ((cfg.get("final_forbidden", {}) or {}).get("motifs", []) or [])
    out: List[str] = []
    out += [str(m).upper() for m in avoid if str(m).strip()]
    out += [str(m).upper() for m in final_forbidden if str(m).strip()]
    return out


def _find_codonopt_output_fasta(out_dir: str) -> str:
    candidates = [
        os.path.join(out_dir, "optimized_sequences.fasta"),
        os.path.join(out_dir, "optimized_sequences.fa"),
        os.path.join(out_dir, "optimized.fasta"),
        os.path.join(out_dir, "results", "optimized_sequences.fasta"),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c

    fastas = sorted(glob.glob(os.path.join(out_dir, "**", "*.fasta"), recursive=True))
    fastas += sorted(glob.glob(os.path.join(out_dir, "**", "*.fa"), recursive=True))
    fastas = [f for f in fastas if os.path.isfile(f)]

    uniq: List[str] = []
    for f in fastas:
        if f not in uniq:
            uniq.append(f)

    if len(uniq) == 1:
        return uniq[0]

    raise FileNotFoundError(
        f"Could not uniquely identify codonopt output FASTA in {out_dir}. Found: {uniq}"
    )


def _inject_codonopt_pythonpath(env: Dict[str, str], workdir: str) -> Dict[str, str]:
    candidates = [
        os.path.join(workdir, "third_party", "codonopt"),
        os.path.join(os.path.dirname(workdir), "third_party", "codonopt"),
    ]
    for c in candidates:
        if os.path.isdir(os.path.join(c, "codonopt")):
            existing = env.get("PYTHONPATH", "")
            env["PYTHONPATH"] = f"{c}{os.pathsep}{existing}" if existing else c
            break
    return env


def run_codonopt(cfg: Dict[str, Any], workdir: str) -> str:
    codonopt_cfg = cfg.get("codonopt", {}) or {}

    out_subdir = codonopt_cfg.get("out_subdir", "codonopt_out")
    out_dir = os.path.join(workdir, out_subdir)
    os.makedirs(out_dir, exist_ok=True)

    seq = codonopt_cfg.get("sequence")
    batch = codonopt_cfg.get("batch_table")
    if bool(seq) == bool(batch):
        raise ValueError("codonopt: provide exactly one of 'sequence' or 'batch_table'.")

    codon_table = codonopt_cfg.get("codon_table")
    if not codon_table:
        raise ValueError("codonopt.codon_table is required.")

    cmd: List[str] = ["python", "-m", "codonopt.main"]

    if seq:
        cmd += ["--sequence", str(seq)]
    if batch:
        cmd += ["--batch-table", str(batch)]

    cmd += ["--codon-table", str(codon_table)]

    if codonopt_cfg.get("codon_table_sheet"):
        cmd += ["--codon-table-sheet", str(codonopt_cfg["codon_table_sheet"])]

    cmd += ["--out", out_dir]

    if codonopt_cfg.get("verbose"):
        cmd += ["--verbose"]
    if codonopt_cfg.get("log_file"):
        cmd += ["--log-file", str(codonopt_cfg["log_file"])]

    for k, flag in [
        ("gc_min", "--gc-min"),
        ("gc_max", "--gc-max"),
        ("max_homopolymer", "--max-homopolymer"),
        ("min_codon_fraction", "--min-codon-fraction"),
        ("seed", "--seed"),
        ("n", "--n"),
    ]:
        if k in codonopt_cfg and codonopt_cfg[k] is not None:
            cmd += [flag, str(codonopt_cfg[k])]

    if codonopt_cfg.get("avoid_motifs"):
        cmd += ["--avoid-motifs", "|".join(map(str, codonopt_cfg["avoid_motifs"]))]

    env = _inject_codonopt_pythonpath(dict(os.environ), workdir)

    proc = subprocess.run(cmd, cwd=workdir, env=env, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "codonopt failed:\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDOUT:\n{proc.stdout}\n"
            f"STDERR:\n{proc.stderr}"
        )

    return _find_codonopt_output_fasta(out_dir)


def _parse_mutation_spec(cfg: Dict[str, Any], baseline_protein: str) -> MutationSpec:
    aa_len = len(baseline_protein)
    mut = cfg.get("mutation")

    if not mut:
        return ("saturation_positions", list(range(1, aa_len + 1)))

    if not isinstance(mut, dict):
        raise ValueError("mutation must be a mapping if provided.")

    if "explicit" in mut and mut["explicit"] is not None:
        return ("explicit", mut["explicit"])

    if "aa_range" in mut and mut["aa_range"] is not None:
        start, end = mut["aa_range"]
        return ("saturation_positions", list(range(int(start), int(end) + 1)))

    if "aa_positions" in mut and mut["aa_positions"] is not None:
        return ("saturation_positions", [int(x) for x in mut["aa_positions"]])

    return ("saturation_positions", list(range(1, aa_len + 1)))


def _normalize_protein(seq: str) -> str:
    s = (seq or "").strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
    if s.endswith("*"):
        s = s[:-1]
    return s


def _base_id(record_id: str) -> str:
    """
    codonopt may emit ids like:
      <original>|job0001|rep001
    Normalize back to <original> for joining with input FASTA IDs.
    """
    rid = (record_id or "").strip()
    if "|job" in rid:
        rid = rid.split("|job", 1)[0]
    if "|rep" in rid:
        rid = rid.split("|rep", 1)[0]
    return rid


def _patch_forced_codons(cds: str, aa_len: int, forced_by_local_aa: Dict[int, str]) -> str:
    """
    Force specific codons at specific AA indices (0-based) in a CDS.
    This is how shuffle guarantees the chosen overhangs are feasible at each junction.
    """
    cds = cds.upper().replace("U", "T").replace(" ", "").replace("\n", "").replace("\r", "")
    if len(cds) != aa_len * 3:
        raise ValueError(f"CDS length mismatch: expected {aa_len*3}, got {len(cds)}")
    arr = list(cds)
    for pos, codon in forced_by_local_aa.items():
        if pos < 0 or pos >= aa_len:
            continue
        i = pos * 3
        arr[i : i + 3] = list(str(codon).upper())
    return "".join(arr)


def run_pipeline(config_path: str, outdir: str) -> None:
    cfg = load_config(config_path)
    os.makedirs(outdir, exist_ok=True)

    workdir = os.path.dirname(os.path.abspath(config_path)) or "."

    inputs = cfg.get("inputs", {}) or {}
    input_fasta = inputs.get("fasta")
    input_type = str(inputs.get("input_type", "protein")).strip().lower()

    if input_type not in ("protein", "cds"):
        raise ValueError("inputs.input_type must be 'protein' or 'cds'.")

    if not input_fasta:
        raise ValueError("inputs.fasta is required.")

    input_fasta_path = os.path.join(workdir, input_fasta)

    codonopt_cfg = cfg.get("codonopt", {}) or {}
    codonopt_enable = bool(codonopt_cfg.get("enable", True))

    gg_cfg = (((cfg.get("golden_gate", {}) or {}).get("fragments", None)) or None)
    if gg_cfg is None:
        raise ValueError("golden_gate.fragments is required.")

    library_mode = ((cfg.get("library_mode", {}) or {}).get("type", "fragment_single_mutation") or "").strip()
    if library_mode != "fragment_single_mutation":
        raise ValueError("This build currently supports library_mode.type = fragment_single_mutation only.")

    # Read original input records for AA-preservation validation
    original_records = read_fasta_records(input_fasta_path)
    original_by_base = {_base_id(rid): seq for (rid, seq) in original_records}

    # Baseline generation (CDS)
    if codonopt_enable:
        if not codonopt_cfg.get("sequence") and not codonopt_cfg.get("batch_table"):
            cfg.setdefault("codonopt", {})
            cfg["codonopt"]["sequence"] = input_fasta

        optimized_fasta = run_codonopt(cfg, workdir)
        optimized = read_fasta_records(optimized_fasta)
    else:
        if input_type != "cds":
            raise ValueError("codonopt.enable=false is only allowed when inputs.input_type='cds'.")
        optimized = read_fasta_records(input_fasta_path)

    # Validate AA preservation for baseline optimized sequences
    for (rid, seq) in optimized:
        base = _base_id(rid)
        if base not in original_by_base:
            raise ValueError(
                f"Record id '{rid}' (base='{base}') is present in codonopt output but not in input FASTA. "
                "IDs must match (or share a common base prefix before '|job'/'|rep') for AA-preservation validation."
            )

        optimized_prot = _normalize_protein(translate_cds(seq))
        if input_type == "protein":
            expected = _normalize_protein(original_by_base[base])
        else:
            expected = _normalize_protein(translate_cds(original_by_base[base]))

        if optimized_prot != expected:
            raise ValueError(
                f"AA sequence mismatch after baseline optimization for record '{rid}' (base='{base}').\n"
                f"Expected (from input): {expected[:120]}{'...' if len(expected) > 120 else ''}\n"
                f"Got (from optimized CDS): {optimized_prot[:120]}{'...' if len(optimized_prot) > 120 else ''}\n"
                "This indicates either an input frame/translation issue or codonopt altered the amino-acid sequence."
            )


    # Codon table used for mutation codon selection + junction feasibility enumeration
    if not codonopt_cfg.get("codon_table"):
        raise ValueError("codonopt.codon_table is required.")

    # NOTE: mutate.py expects a mapping AA -> list of codon entries (CodonEntry is supported if you patched mutate.py)
    table = load_codon_usage_xlsx(
        codonopt_cfg["codon_table"],
        codonopt_cfg.get("codon_table_sheet"),
    )

    forbidden = _collect_forbidden_motifs(cfg)

    # Outputs
    out_cfg = cfg.get("outputs", {}) or {}
    baseline_out = os.path.join(outdir, out_cfg.get("baseline_cds_fasta", "baseline_optimized.fasta"))
    variants_out = os.path.join(outdir, out_cfg.get("variants_cds_fasta", "variants.fasta"))

    wt_frags_path = os.path.join(outdir, out_cfg.get("wt_fragments_fasta", "wt_fragments_with_gg.fasta"))
    mut_frags_path = os.path.join(outdir, out_cfg.get("mut_fragments_fasta", "mut_fragments_with_gg.fasta"))
    report_path = os.path.join(outdir, out_cfg.get("report_tsv", "report.tsv"))
    assembly_plan_path = os.path.join(outdir, out_cfg.get("assembly_plan_tsv", "assembly_plan.tsv"))
    pools_dir = os.path.join(outdir, out_cfg.get("pools_dir", "pools"))
    os.makedirs(pools_dir, exist_ok=True)

    # Collections
    variant_entries: List[Tuple[str, str]] = []
    wt_fragment_entries: List[Tuple[str, str]] = []
    mut_fragment_entries: List[Tuple[str, str]] = []

    pools: Dict[int, List[Tuple[str, str]]] = {}
    report_rows: List[Dict[str, Any]] = []
    assembly_rows: List[Dict[str, Any]] = []

    # Persist baseline optimized CDS
    write_fasta_records(baseline_out, [(rid, seq) for (rid, seq) in optimized])

    for (gene_id_full, baseline_cds) in optimized:
        gene_id = _base_id(gene_id_full)

        baseline_prot = translate_cds(baseline_cds)

        mut_mode, mut_payload = _parse_mutation_spec(cfg, baseline_prot)
        cut_points, enforce_frame = _merge_cut_points(cfg, gene_id)

        # Build core fragments (before forcing junction codons)
        wt_frags_core = split_cds_to_fragment_dicts(
            baseline_cds,
            aa_cut_points=cut_points,
            enforce_frame=enforce_frame,
        )
        n_frags = len(wt_frags_core)
        if n_frags < 1:
            raise RuntimeError(f"{gene_id}: no fragments produced. Check cut points.")

        # --- Golden Gate packaging (converge on SHUFFLE logic) ---
        enzyme = str(gg_cfg.get("enzyme") or "BsaI")
        spacer_left = str(gg_cfg.get("spacer_left") or "")
        spacer_right = str(gg_cfg.get("spacer_right") or "")

        aa_records = [(gene_id, baseline_prot)]
        cut_points_by_gene = {gene_id: cut_points}

        # Choose overhangs globally using shuffle feasibility + matrix crosstalk minimization
        overhang_by_pos, forced_codons_by_gene = select_global_overhangs_and_forced_codons(
            cfg=cfg,
            aa_records=aa_records,
            n_fragments=n_frags,
            cut_points_by_gene=cut_points_by_gene,
        )

        # ---- Overhang audit output + crosstalk threshold ----
        from frogger.overhangs import crosstalk_rows
        from frogger.common.io import write_tsv

        gg = cfg.get("golden_gate", {}) or {}
        assign = gg.get("_selected_overhangs_by_pos", overhang_by_pos)
        worst = gg.get("_selected_overhangs_worst_score", None)

        selected_rows = [{"pos_index": int(pos), "overhang": str(oh)} for pos, oh in sorted(assign.items())]
        write_tsv(Path(outdir) / "overhangs_selected.tsv", selected_rows)

        df_matrix = load_pairing_matrix_xlsx(
            Path(gg["overhang_matrix_xlsx"]),
            sheet_name=gg.get("overhang_matrix_sheet"),
        )
        write_tsv(Path(outdir) / "overhang_crosstalk.tsv", crosstalk_rows(df_matrix, assign))

        max_allowed = float(gg.get("max_worst_crosstalk", 0.1))
        if worst is not None and float(worst) > max_allowed:
            raise RuntimeError(
                f"Overhang assignment worst cross-talk {float(worst):.6g} exceeds max_worst_crosstalk={max_allowed:.6g}. "
                f"See {Path(outdir) / 'overhang_crosstalk.tsv'}"
            )

        forced_map = forced_codons_by_gene.get(gene_id, {})
        baseline_cds_for_gg = _patch_forced_codons(
            baseline_cds,
            aa_len=len(baseline_prot),
            forced_by_local_aa=forced_map,
        )

        # Re-split on the patched CDS so cores reflect forced codons
        wt_frags_core = split_cds_to_fragment_dicts(
            baseline_cds_for_gg,
            aa_cut_points=cut_points,
            enforce_frame=enforce_frame,
        )

        overhangs_by_pos = compute_overhangs_by_frag_pos(overhang_by_pos, n_fragments=n_frags)
        primers_by_pos = build_primers_by_pos_for_fragments(cfg, fragments_cds=wt_frags_core, n_fragments=n_frags)

        wt_frags_clone = build_cloning_sequences_for_fragments(
            fragments=wt_frags_core,
            enzyme=enzyme,
            overhangs_by_pos=overhangs_by_pos,
            primers_by_pos=primers_by_pos,
            spacer_left=spacer_left,
            spacer_right=spacer_right,
        )

        # WT fragments + seed pools with WT
        for f in wt_frags_clone:
            wt_id = f"{gene_id}|WT|frag{int(f['frag_index']):02d}"
            wt_fragment_entries.append((wt_id, str(f["cloning_seq"])))

        for i in range(1, n_frags + 1):
            pools.setdefault(i, [])
            for f in wt_frags_clone:
                pools[i].append((f"{gene_id}|WT|frag{int(f['frag_index']):02d}", str(f["cloning_seq"])))

        # Variants iterator
        if mut_mode == "explicit":
            variants_iter = generate_explicit_variants(
                gene_id=gene_id,
                baseline_cds=baseline_cds_for_gg,  # baseline already forced at junctions
                forbidden_motifs=forbidden,
                codon_table=table,
                explicit_mutations=mut_payload,  # type: ignore[arg-type]
            )
        else:
            variants_iter = generate_saturation_variants(
                gene_id=gene_id,
                baseline_cds=baseline_cds_for_gg,  # baseline already forced at junctions
                forbidden_motifs=forbidden,
                codon_table=table,
                include_wt=True,
                positions_to_mutate=mut_payload,  # type: ignore[arg-type]
            )

        # Map AA position -> fragment index using aa_start0/aa_end0
        for var in variants_iter:
            assert isinstance(var, Variant)

            variant_id = f"{gene_id}|{var.mutation_label}"
            variant_entries.append((variant_id, var.cds))

            mut_aa0 = var.pos1 - 1
            mutated_frag_index: Optional[int] = None
            for core in wt_frags_core:
                if int(core["aa_start0"]) <= mut_aa0 < int(core["aa_end0"]):
                    mutated_frag_index = int(core["frag_index"])
                    break

            if mutated_frag_index is None:
                raise RuntimeError(
                    f"{gene_id}|{var.mutation_label}: failed to map mutation aa0={mut_aa0} to a fragment "
                    f"(n_frags={len(wt_frags_core)})."
                )

            # Patch variant CDS at junction codons to keep chosen overhangs feasible
            var_cds_for_gg = _patch_forced_codons(
                var.cds,
                aa_len=len(baseline_prot),
                forced_by_local_aa=forced_map,
            )

            var_frags_core = split_cds_to_fragment_dicts(
                var_cds_for_gg,
                aa_cut_points=cut_points,
                enforce_frame=enforce_frame,
            )
            var_frags_clone = build_cloning_sequences_for_fragments(
                fragments=var_frags_core,
                enzyme=enzyme,
                overhangs_by_pos=overhangs_by_pos,
                primers_by_pos=primers_by_pos,
                spacer_left=spacer_left,
                spacer_right=spacer_right,
            )

            mut_frag = next(f for f in var_frags_clone if int(f["frag_index"]) == mutated_frag_index)
            mut_frag_id = f"{gene_id}|{var.mutation_label}|frag{mutated_frag_index:02d}"

            mut_fragment_entries.append((mut_frag_id, str(mut_frag["cloning_seq"])))
            pools[mutated_frag_index].append((mut_frag_id, str(mut_frag["cloning_seq"])))

            report_rows.append(
                {
                    "gene_id": gene_id,
                    "variant_id": variant_id,
                    "mutation": var.mutation_label,
                    "position_aa": var.pos1,
                    "wt_aa": var.wt_aa,
                    "target_aa": var.target_aa,
                    "codon_used": var.codon_used,
                    "mutated_frag_index": mutated_frag_index,
                    "n_frags": n_frags,
                }
            )

            row: Dict[str, Any] = {
                "gene_id": gene_id,
                "variant_id": variant_id,
                "mutation": var.mutation_label,
                "mutated_frag_index": mutated_frag_index,
            }
            for f in wt_frags_clone:
                frag_index = int(f["frag_index"])
                frag_key = f"frag{frag_index:02d}_id"
                if frag_index == mutated_frag_index:
                    row[frag_key] = mut_frag_id
                else:
                    row[frag_key] = f"{gene_id}|WT|frag{frag_index:02d}"
            assembly_rows.append(row)

    # Write FASTA outputs
    write_fasta_records(variants_out, variant_entries)
    write_fasta_records(wt_frags_path, wt_fragment_entries)
    write_fasta_records(mut_frags_path, mut_fragment_entries)

    # Write pools
    for idx, entries in sorted(pools.items(), key=lambda kv: kv[0]):
        seen = set()
        uniq = []
        for h, s in entries:
            key = (h, s)
            if key not in seen:
                uniq.append((h, s))
                seen.add(key)
        write_fasta_records(os.path.join(pools_dir, f"pool_frag{idx:02d}.fasta"), uniq)

    # Write report.tsv
    report_header = [
        "gene_id",
        "variant_id",
        "mutation",
        "position_aa",
        "wt_aa",
        "target_aa",
        "codon_used",
        "mutated_frag_index",
        "n_frags",
    ]
    write_tsv(report_path, report_rows)

    # Write assembly_plan.tsv
    max_idx = 0
    for r in assembly_rows:
        for k in r.keys():
            if k.startswith("frag") and k.endswith("_id"):
                try:
                    num = int(k[4:6])
                    max_idx = max(max_idx, num)
                except Exception:
                    pass

    plan_header = ["gene_id", "variant_id", "mutation", "mutated_frag_index"] + [
        f"frag{i:02d}_id" for i in range(1, max_idx + 1)
    ]
    write_tsv(assembly_plan_path, assembly_rows)

def run(config_path: str, outdir: str) -> None:
    """
    CLI entrypoint shim. Keep this name stable across modules.
    """
    return run_pipeline(config_path=config_path, outdir=outdir)
