from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from Bio.Seq import Seq

from frogger.util import effective_seed, run_codonopt_for_fragment
from frogger.common.gg.splitters import split_aa_by_cut_points
from frogger.common.codon_table import load_codon_usage_xlsx
from frogger.junction_overhangs import enumerate_boundary_overhangs
from frogger.overhangs import load_pairing_matrix_xlsx, choose_overhangs_by_position, OverhangFilters
from frogger.gg_adapters import build_cloning_sequences_for_fragments
from frogger.recombine import generate_constructs
from frogger.checks import check_construct
from frogger.report import build_report_row, build_hit_rows
from frogger.common.io import read_fasta_records, write_fasta_records, write_tsv



def _load_aa_inputs(cfg: dict) -> List[tuple[str, str]]:
    inp = cfg.get("inputs", {}) or {}
    aa_path = inp.get("aa_fasta") or inp.get("protein_fasta") or inp.get("fasta")
    if not aa_path:
        aa_path = (cfg.get("codonopt", {}) or {}).get("sequence")
    if not aa_path:
        raise ValueError("No AA FASTA provided. Set inputs.aa_fasta (preferred).")
    return read_fasta_records(Path(aa_path))


def _get_cut_points_aa(cfg: dict, gene_id: str) -> List[int]:
    split_cfg = cfg["splits"]
    global_cuts = split_cfg.get("global_cut_points")
    per_gene = split_cfg.get("per_gene_cut_points", {}) or {}
    cut_points = per_gene.get(gene_id, global_cuts)
    if cut_points is None:
        raise ValueError("No cut points set (splits.global_cut_points or splits.per_gene_cut_points).")
    return [int(x) for x in list(cut_points)]


def _patch_forced_codons(cds: str, aa_len: int, forced_by_local_aa: Dict[int, str]) -> str:
    cds = cds.upper().replace("U", "T")
    if len(cds) != aa_len * 3:
        raise ValueError(f"CDS length mismatch: expected {aa_len*3}, got {len(cds)}")
    arr = list(cds)
    for pos, codon in forced_by_local_aa.items():
        if pos < 0 or pos >= aa_len:
            continue
        i = pos * 3
        arr[i:i+3] = list(codon.upper())
    return "".join(arr)


def _select_global_overhangs_and_forced_codons(
    cfg: dict,
    aa_records: List[Tuple[str, str]],
    n_fragments: int,
    cut_points_by_gene: Dict[str, List[int]],
) -> Tuple[Dict[int, str], Dict[str, Dict[int, str]]]:
    """
    Returns:
      - overhang_by_pos: {pos_index(0..k): 'ACGT'}
      - forced_codons_by_gene: {gene_id: {aa_index: 'ATG', ...}}
    """
    gg = cfg.get("golden_gate", {}) or {}
    matrix_path = gg.get("overhang_matrix_xlsx")
    if not matrix_path:
        raise ValueError("golden_gate.overhang_matrix_xlsx is required.")
    df = load_pairing_matrix_xlsx(Path(matrix_path), sheet_name=gg.get("overhang_matrix_sheet"))

    # Codon usage table (for feasible overhang enumeration)
    c_cfg = cfg.get("codonopt", {}) or {}
    codon_table_path = c_cfg.get("codon_table")
    if not codon_table_path:
        raise ValueError("codonopt.codon_table is required (used to enumerate synonymous codons near junctions).")
    codon_map = load_codon_usage_xlsx(Path(codon_table_path), sheet_name=c_cfg.get("codon_table_sheet"))
    avoid_codons = c_cfg.get("avoid_codons", []) or []

    # For each junction position, compute candidate overhangs that are feasible for ALL genes
    # Positions: 0 = N-term, 1..k-1 internal boundaries, k = C-term
    candidates_by_pos: Dict[int, List[str]] = {}
    per_gene_choices: Dict[str, Dict[int, Dict[str, object]]] = {}

    for gene_id, aa in aa_records:
        aa = aa.strip()
        cuts = cut_points_by_gene[gene_id]
        if len(aa) < 2:
            raise ValueError(f"Sequence '{gene_id}' too short for terminus overhang search (<2 aa).")

        boundaries: List[Tuple[int, str, str]] = []
        # N-term uses aa[0], aa[1]
        boundaries.append((0, aa[0], aa[1]))
        # internal boundaries
        for idx, cp in enumerate(cuts, start=1):
            if cp <= 0 or cp >= len(aa):
                raise ValueError(f"Cut point {cp} out of bounds for '{gene_id}' length {len(aa)}")
            boundaries.append((idx, aa[cp - 1], aa[cp]))
        # C-term uses last 2 aa
        boundaries.append((len(cuts) + 1, aa[-2], aa[-1]))

        per_gene_choices[gene_id] = {}
        for pos_index, aL, aR in boundaries:
            cand_map = enumerate_boundary_overhangs(aL, aR, codon_map, avoid_codons=avoid_codons)
            per_gene_choices[gene_id][pos_index] = cand_map

    # Intersect candidates across genes for each position
    for pos_index in range(0, n_fragments + 1):
        keys = None
        for gene_id, _aa in aa_records:
            gene_cands = set(per_gene_choices[gene_id][pos_index].keys())
            keys = gene_cands if keys is None else (keys & gene_cands)
        if not keys:
            raise RuntimeError(
                f"No feasible 4-mer overhangs exist at junction position {pos_index} "
                f"for all genes given the fixed cut points and codon table."
            )
        candidates_by_pos[pos_index] = sorted(keys)

    # Choose one overhang per position minimizing cross-talk
    flt_cfg = (gg.get("overhang_filters", {}) or {})
    flt = OverhangFilters(
        disallow_palindromes=bool(flt_cfg.get("disallow_palindromes", True)),
        disallow_revcomp_self=bool(flt_cfg.get("disallow_revcomp_self", True)),
        min_gc=float(flt_cfg.get("min_gc", 0.25)),
        max_gc=float(flt_cfg.get("max_gc", 0.75)),
        max_homopolymer=int(flt_cfg.get("max_homopolymer", 3)),
    )

    assign, worst = choose_overhangs_by_position(
        df=df,
        candidates_by_pos=candidates_by_pos,
        filters=flt,
        beam_width=int(gg.get("beam_width", 200)),
    )

    gg["_selected_overhangs_by_pos"] = assign
    gg["_selected_overhangs_worst_score"] = worst
    cfg["golden_gate"] = gg

    # Now, for each gene, choose concrete codon pairs that realize the selected overhang at each position,
    # preferring higher codon usage product.
    forced_codons_by_gene: Dict[str, Dict[int, str]] = {}
    for gene_id, aa in aa_records:
        aa = aa.strip()
        cuts = cut_points_by_gene[gene_id]
        forced: Dict[int, str] = {}

        # Helper to map position index -> aa indices involved
        def aa_indices_for_pos(pos_index: int) -> Tuple[int, int]:
            if pos_index == 0:
                return 0, 1
            if pos_index == n_fragments:
                return len(aa) - 2, len(aa) - 1
            # internal pos: pos_index corresponds to cut point at cuts[pos_index-1]
            cp = cuts[pos_index - 1]
            return cp - 1, cp

        for pos_index in range(0, n_fragments + 1):
            oh = assign[pos_index]
            choices = per_gene_choices[gene_id][pos_index]
            if oh not in choices:
                raise RuntimeError(f"Internal error: chosen overhang {oh} not feasible for gene {gene_id} pos {pos_index}")
            choice = choices[oh]
            iL, iR = aa_indices_for_pos(pos_index)
            # Force codons for both amino acids at this boundary
            forced[iL] = choice.left_codon
            forced[iR] = choice.right_codon

        forced_codons_by_gene[gene_id] = forced

    return assign, forced_codons_by_gene


def run_pipeline(cfg: dict, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)

    seed_codonopt = effective_seed(cfg["codonopt"].get("seed"), cfg["defaults"]["default_seed"])
    seed_reasm = effective_seed(cfg["reassembly"].get("seed"), cfg["defaults"]["default_seed"])

    aa_records = _load_aa_inputs(cfg)

    # Precompute AA fragments and cut points per gene
    fragments_aa: List[Dict] = []
    cut_points_by_gene: Dict[str, List[int]] = {}
    for gene_id, aa in aa_records:
        cut_points = _get_cut_points_aa(cfg, gene_id)
        cut_points_by_gene[gene_id] = cut_points
        frags = split_aa_by_cut_points(aa, cut_points)
        for i, (a, b, frag_aa) in enumerate(frags, start=1):
            fragments_aa.append(
                {
                    "gene_id": gene_id,
                    "frag_index": i,
                    "aa_start": a,
                    "aa_end": b,
                    "aa_seq": frag_aa,
                }
            )
    n_fragments = max(int(f["frag_index"]) for f in fragments_aa)

    # Choose GLOBAL overhangs per junction position (0..k) AND forced codons per gene/AA-index
    overhang_by_pos, forced_codons_by_gene = _select_global_overhangs_and_forced_codons(
        cfg=cfg,
        aa_records=aa_records,
        n_fragments=n_fragments,
        cut_points_by_gene=cut_points_by_gene,
    )

    # ---- Overhang audit output + crosstalk threshold ----
    from frogger.overhangs import crosstalk_rows  # uses the pairing matrix
    from frogger.common.io import write_tsv

    gg = cfg.get("golden_gate", {}) or {}
    # Prefer the stored selection (set inside the selector), but fall back to the returned mapping.
    assign = gg.get("_selected_overhangs_by_pos", overhang_by_pos)
    worst = gg.get("_selected_overhangs_worst_score", None)

    # Write selected overhangs (pos 0..n_fragments)
    selected_rows = [{"pos_index": int(pos), "overhang": str(oh)} for pos, oh in sorted(assign.items())]
    write_tsv(Path(outdir) / "overhangs_selected.tsv", selected_rows)

    # Write full pairwise cross-talk table
    df_matrix = load_pairing_matrix_xlsx(
        Path(gg["overhang_matrix_xlsx"]),
        sheet_name=gg.get("overhang_matrix_sheet"),
    )
    write_tsv(Path(outdir) / "overhang_crosstalk.tsv", crosstalk_rows(df_matrix, assign))

    # Enforce threshold (default 0.1)
    max_allowed = float(gg.get("max_worst_crosstalk", 0.1))
    if worst is not None and float(worst) > max_allowed:
        raise RuntimeError(
            f"Overhang assignment worst cross-talk {float(worst):.6g} exceeds max_worst_crosstalk={max_allowed:.6g}. "
            f"See {Path(outdir) / 'overhang_crosstalk.tsv'}"
        )

    # Map junction overhangs into fragment-end overhangs
    overhangs_by_frag_pos: Dict[int, Dict[str, str]] = {}
    for i in range(1, n_fragments + 1):
        overhangs_by_frag_pos[i] = {"left": overhang_by_pos[i - 1], "right": overhang_by_pos[i]}

    # Codonopt per fragment, then patch forced codons for AA positions that fall inside this fragment
    fragments_cds: List[Dict] = []
    c_cfg = cfg["codonopt"]
    for f in fragments_aa:
        gene_id = f["gene_id"]
        frag_id = f'{gene_id}__pos{int(f["frag_index"]):02d}__aa{f["aa_start"]}-{f["aa_end"]}'
        optimized_cds, _, _ = run_codonopt_for_fragment(
            codonopt_cfg=c_cfg,
            outdir=outdir,
            seed=int(seed_codonopt),
            fragment_id=frag_id,
            aa_seq=f["aa_seq"],
        )

        # Patch codons inside this fragment that are forced by junction overhang selection
        forced_global = forced_codons_by_gene.get(gene_id, {})
        forced_local: Dict[int, str] = {}
        for aa_idx, codon in forced_global.items():
            if int(f["aa_start"]) <= aa_idx < int(f["aa_end"]):
                forced_local[aa_idx - int(f["aa_start"])] = codon

        patched_cds = _patch_forced_codons(optimized_cds, aa_len=len(f["aa_seq"]), forced_by_local_aa=forced_local)

        # sanity translation
        if str(Seq(patched_cds).translate(to_stop=False)) != f["aa_seq"]:
            raise RuntimeError(
                f"Forced-codon patching changed translation for {frag_id}. "
                f"This indicates an internal bug in codon forcing."
            )

        fragments_cds.append(
            {
                "fragment_id": frag_id,                 # NEW (stable ID for reporting)
                "gene_id": gene_id,
                "frag_index": f["frag_index"],
                "start": f["aa_start"],
                "end": f["aa_end"],
                "core_seq": patched_cds,
                "aa_seq": f["aa_seq"],
                "forced_codons_local": forced_local,    # NEW (needed for repair)
                "aa_len": len(f["aa_seq"]),             # NEW (optional but handy)
            }
        )
    
    # --- QC + repair pass on fragment core_seq (optional) -----------------------
    repair_cfg = (cfg.get("repair", {}) or {})
    if bool(repair_cfg.get("enable", False)):
        from frogger.common.repair import repair_fragments_loop

        qc_cfg = (repair_cfg.get("qc", {}) or {})
        max_rounds = int(repair_cfg.get("max_rounds", 1))

        def _translate_nt(nt: str) -> str:
            aa = str(Seq(nt).translate(to_stop=False))
            # fragments should not contain stops; if they do, keep '*' so mismatch is caught
            return aa

        def _codonopt_runner(fragment_dict: dict):
            # Re-run codonopt on the AA fragment only (deterministic seed)
            new_nt, _, _ = run_codonopt_for_fragment(
                codonopt_cfg=c_cfg,
                outdir=outdir,
                seed=int(seed_codonopt),
                fragment_id=str(fragment_dict["fragment_id"]) + f"__repair{1}",
                aa_seq=str(fragment_dict["aa_seq"]),
            )
            # return (new_nt, codonopt_version)
            return new_nt, ""  # if you have version, return it here

        fragments_cds, repair_rows = repair_fragments_loop(
            fragments_cds,
            qc_cfg=qc_cfg,
            enable=True,
            max_rounds=max_rounds,
            codonopt_fragment_runner=_codonopt_runner,
            aa_translate=_translate_nt,
            patch_forced_codons=_patch_forced_codons,
        )

        if repair_rows:
            repair_path = outdir / cfg["outputs"].get("repair_report_tsv", "repair_report.tsv")
            write_tsv(repair_path, repair_rows)
    # --------------------------------------------------------------------------

    from frogger.primers import build_avoid_kmers, generate_primers_by_pos

    primers_cfg = cfg.get("primers", {}) or {}
    mode = (primers_cfg.get("mode") or "fixed").lower()

    if mode == "fixed":
        primers_by_pos = primers_cfg.get("fixed_by_pos") or {}
    elif mode == "generate":
        g = primers_cfg.get("generate", {}) or {}
        k = int(g.get("max_internal_match", 12))
        avoid = build_avoid_kmers([f["core_seq"] for f in fragments_cds], k=k)
        primers_by_pos = generate_primers_by_pos(
            n_positions=n_fragments,
            length=int(g.get("length", 18)),
            avoid_kmers=avoid,
            k=k,
            min_gc=float(g.get("min_gc", 0.35)),
            max_gc=float(g.get("max_gc", 0.65)),
            max_homopolymer=int(g.get("max_homopolymer", 3)),
            avoid_motifs=list(g.get("avoid_motifs", [])),
            seed=int(cfg.get("defaults", {}).get("default_seed", 1337)),
        )
    else:
        raise ValueError(f"Unknown primers.mode: {mode}")

    # Build cloning sequences (primer + TypeIIS + OH + core + OH + TypeIIS + primer)
    gg = cfg.get("golden_gate", {}) or {}
    enzyme = gg.get("enzyme", "BsaI")
    spacer_left = gg.get("spacer_left", "")
    spacer_right = gg.get("spacer_right", "")

    gg_frags = build_cloning_sequences_for_fragments(
        fragments=fragments_cds,
        enzyme=enzyme,
        overhangs_by_pos=overhangs_by_frag_pos,
        primers_by_pos=primers_by_pos,
        spacer_left=spacer_left,
        spacer_right=spacer_right,
    )

    # Outputs: fragments FASTA + fragments table
    fragments_fasta_path = outdir / cfg["outputs"]["fragments_fasta"]
    frag_records = []
    table_rows = []
    for f in gg_frags:
        header = f'{f["gene_id"]}|pos={f["frag_index"]}|aa={f["start"]}-{f["end"]}'
        frag_records.append((header, f["cloning_seq"]))
        table_rows.append(
            {
                "gene_id": f["gene_id"],
                "frag_index": int(f["frag_index"]),
                "aa_start": int(f["start"]),
                "aa_end": int(f["end"]),
                "enzyme": enzyme,
                "left_overhang": overhangs_by_frag_pos[int(f["frag_index"])]["left"],
                "right_overhang": overhangs_by_frag_pos[int(f["frag_index"])]["right"],
                "left_adapter": f["left_adapter"],
                "right_adapter": f["right_adapter"],
                "left_primer": f["primer_left"],
                "right_primer": f["primer_right"],
                "core_cds": f["core_seq"],
                "cloning_seq": f["cloning_seq"],
            }
        )
    write_fasta_records(fragments_fasta_path, frag_records)

    table_path = outdir / cfg["outputs"].get("fragments_table_csv", "fragments_table.csv")
    pd.DataFrame(table_rows).to_csv(table_path, index=False)

    # Recombine fragments across genes (existing behavior)
    re_cfg = cfg["reassembly"]
    constructs = generate_constructs(
        fragments=fragments_cds,
        max_constructs=int(re_cfg["max_constructs"]),
        sample_n=int(re_cfg["sample_n"]),
        seed=int(seed_reasm),
        deduplicate=bool(re_cfg.get("deduplicate", True)),
    )

    final_cfg = cfg.get("final_checks", {}) or {}
    extra_final_motifs = (cfg.get("final_forbidden", {}) or {}).get("motifs", []) or []

    require_no_internal_stops = bool(final_cfg.get("require_no_internal_stops", True))
    require_len_multiple_of_3 = bool(final_cfg.get("require_length_multiple_of_3", True))

    report_rows = []
    hit_rows = []
    orf_records = []

    codonopt_version = "v1.2"
    for c in constructs:
        result = check_construct(
            construct=c,
            codonopt_cfg=cfg["codonopt"],
            final_forbidden_motifs=extra_final_motifs,
            require_no_internal_stops=require_no_internal_stops,
            require_len_multiple_of_3=require_len_multiple_of_3,
        )

        report_rows.append(build_report_row(c, result, codonopt_version, seed_codonopt, seed_reasm))
        hit_rows.extend(build_hit_rows(c, result))

        status = "PASS" if result["passes_all"] else "FAIL"
        orf_records.append((f'{c["construct_id"]}|status={status}|sources={c["frag_sources_str"]}', c["core_orf_seq"]))

    reassembled_fasta_path = outdir / cfg["outputs"]["reassembled_fasta"]
    write_fasta_records(reassembled_fasta_path, orf_records)

    report_path = outdir / cfg["outputs"]["report_tsv"]
    write_tsv(report_path, report_rows)

    hits_path = outdir / cfg["outputs"]["junction_hits_tsv"]
    write_tsv(hits_path, hit_rows)
