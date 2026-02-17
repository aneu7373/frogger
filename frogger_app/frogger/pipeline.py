from __future__ import annotations

from collections import Counter
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
        arr[i : i + 3] = list(codon.upper())
    return "".join(arr)


def _require_terminal_overhangs(cfg: dict) -> Tuple[str, str]:
    """
    (2) Terminal overhangs are REQUIRED. If not provided, raise a hard error.
    Supports legacy key names too.
    """
    gg = cfg.get("golden_gate", {}) or {}
    fr = gg.get("fragments", {}) or {}

    oh5 = fr.get("terminal_overhang_5p", fr.get("term_overhang_5p"))
    oh3 = fr.get("terminal_overhang_3p", fr.get("term_overhang_3p"))

    if not oh5 or not oh3:
        raise ValueError(
            "Terminal overhangs are required. Please set "
            "golden_gate.fragments.terminal_overhang_5p and terminal_overhang_3p "
            "(or legacy keys term_overhang_5p/term_overhang_3p)."
        )
    oh5 = str(oh5).upper().strip()
    oh3 = str(oh3).upper().strip()
    if len(oh5) != 4 or len(oh3) != 4:
        raise ValueError("Terminal overhangs must be 4 bp each.")
    return oh5, oh3


def _apply_consensus_downstream_aas(
    aa_records: List[Tuple[str, str]],
    cut_points_by_gene: Dict[str, List[int]],
) -> List[Tuple[str, str]]:
    """
    (5) If amino acids differ at a cut site between genes, use the most common in the set.
        OK if this changes sequences (shuffle context).

    Clarification with (1):
      - We search the FIRST TWO amino acids AFTER the cut site.
      - For a cut point cp (boundary between cp-1 | cp), downstream AAs are:
          aa1 = aa[cp]
          aa2 = aa[cp+1]
      - We compute consensus aa1/aa2 across genes and patch each AA record accordingly.
    """
    if len(aa_records) <= 1:
        return aa_records

    # Use the first gene as reference for number of cut points.
    ref_gene = aa_records[0][0]
    n_cuts = len(cut_points_by_gene[ref_gene])

    # Compute consensus for each cut index.
    consensus: Dict[int, Tuple[str, str]] = {}
    for cut_idx in range(n_cuts):
        c1 = Counter()
        c2 = Counter()
        for gene_id, aa in aa_records:
            aa = aa.strip()
            cuts = cut_points_by_gene[gene_id]
            cp = int(cuts[cut_idx])
            if cp < 0 or cp + 1 >= len(aa):
                raise ValueError(
                    f"Cut point {cp} invalid for downstream-AA window in '{gene_id}' length {len(aa)}. "
                    "Cut points must allow two amino acids after the cut (cp <= len(aa)-2)."
                )
            c1[aa[cp]] += 1
            c2[aa[cp + 1]] += 1

        # deterministic tiebreaker
        aa1 = sorted(c1.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        aa2 = sorted(c2.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        consensus[cut_idx] = (aa1, aa2)

    # Apply patches
    patched: List[Tuple[str, str]] = []
    for gene_id, aa in aa_records:
        aa = aa.strip()
        arr = list(aa)
        cuts = cut_points_by_gene[gene_id]
        for cut_idx, (aa1, aa2) in consensus.items():
            cp = int(cuts[cut_idx])
            arr[cp] = aa1
            arr[cp + 1] = aa2
        patched.append((gene_id, "".join(arr)))
    return patched


def _select_global_overhangs_and_forced_codons(
    cfg: dict,
    aa_records: List[Tuple[str, str]],
    n_fragments: int,
    cut_points_by_gene: Dict[str, List[int]],
) -> Tuple[Dict[int, str], Dict[str, Dict[int, str]], List[Tuple[str, str]]]:
    """
    Returns:
      - overhang_by_pos: {pos_index(0..k): 'ACGT'}
      - forced_codons_by_gene: {gene_id: {aa_index: 'ATG', ...}}  (internal junctions only)
      - aa_records_patched: AA records possibly modified by consensus junction AA patching (shuffle)

    Implements requested updates:
      (1) Search FIRST TWO amino acids AFTER cut site (downstream) for internal junction overhang feasibility.
      (2) 5' and 3' terminal overhangs are REQUIRED and hard-coded (not chosen by search).
      (3) Compatibility scoring includes reverse complements (handled in frogger/overhangs.py).
          Also double-check 3' block uses revcomp(overhang) (handled in gg_adapters.build_cloning_seq).
      (4) Restrict codon search to top-K frequent codons per AA (default K=2; pulled from codon table entries).
      (5) If amino acids differ at selected site between genes, use most common in the set (patch AA records).
    """
    # Required terminal overhangs
    term5, term3 = _require_terminal_overhangs(cfg)

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

    # (4) top-k codons per AA for enumeration (default k=2)
    search_cfg = (gg.get("overhang_search", {}) or {})
    top_k = int(search_cfg.get("top_codons_per_aa", 2))
    if top_k < 1:
        top_k = 1

    # (5) Patch AA records at junction downstream AAs to consensus across genes (shuffle assumption)
    aa_records_patched = _apply_consensus_downstream_aas(aa_records, cut_points_by_gene)

    # For each INTERNAL junction position (1..k-1), compute candidate overhangs feasible for ALL genes
    # Positions:
    #   0 = N-term fixed (term5)
    #   1..k-1 internal boundaries chosen by search
    #   k = C-term fixed (term3)
    candidates_by_pos: Dict[int, List[str]] = {}
    per_gene_choices: Dict[str, Dict[int, Dict[str, object]]] = {}

    for gene_id, aa in aa_records_patched:
        aa = aa.strip()
        cuts = cut_points_by_gene[gene_id]
        if len(cuts) != (n_fragments - 1):
            raise ValueError(
                f"Gene '{gene_id}' has {len(cuts)} cut points but n_fragments={n_fragments} "
                f"(expected {n_fragments-1} cut points)."
            )

        # validate that downstream two-AA window exists for every cut
        for cp in cuts:
            if cp < 0 or cp + 1 >= len(aa):
                raise ValueError(
                    f"Cut point {cp} invalid for downstream-AA window in '{gene_id}' length {len(aa)}. "
                    "Cut points must allow two amino acids after the cut (cp <= len(aa)-2)."
                )

        per_gene_choices[gene_id] = {}
        for pos_index, cp in enumerate(cuts, start=1):
            # (1) downstream two AAs after cut: aa[cp], aa[cp+1]
            a1 = aa[cp]
            a2 = aa[cp + 1]
            cand_map = enumerate_boundary_overhangs(
                a1,
                a2,
                codon_map,
                avoid_codons=avoid_codons,
                top_k_codons=top_k,
            )
            per_gene_choices[gene_id][pos_index] = cand_map

    # Intersect candidates across genes for each INTERNAL position
    for pos_index in range(1, n_fragments):
        keys = None
        for gene_id, _aa in aa_records_patched:
            gene_cands = set(per_gene_choices[gene_id][pos_index].keys())
            keys = gene_cands if keys is None else (keys & gene_cands)
        if not keys:
            raise RuntimeError(
                f"No feasible 4-mer overhangs exist at internal junction position {pos_index} "
                f"for all genes given the fixed cut points and codon table."
            )
        candidates_by_pos[pos_index] = sorted(keys)

    # Choose one overhang per INTERNAL position minimizing cross-talk (RC-aware scoring in overhangs.py)
    flt_cfg = (gg.get("overhang_filters", {}) or {})
    flt = OverhangFilters(
        disallow_palindromes=bool(flt_cfg.get("disallow_palindromes", True)),
        disallow_revcomp_self=bool(flt_cfg.get("disallow_revcomp_self", True)),
        min_gc=float(flt_cfg.get("min_gc", 0.25)),
        max_gc=float(flt_cfg.get("max_gc", 0.75)),
        max_homopolymer=int(flt_cfg.get("max_homopolymer", 3)),
    )

    internal_assign, worst = choose_overhangs_by_position(
        df=df,
        candidates_by_pos=candidates_by_pos,
        filters=flt,
        beam_width=int(gg.get("beam_width", 200)),
    )

    # Full assignment includes required terminal overhangs
    assign: Dict[int, str] = {0: term5, n_fragments: term3}
    for pos_index, oh in internal_assign.items():
        assign[int(pos_index)] = str(oh).upper()

    gg["_selected_overhangs_by_pos"] = assign
    gg["_selected_overhangs_worst_score"] = worst
    cfg["golden_gate"] = gg

    # For each gene, choose concrete codon pairs that realize the selected overhang at each INTERNAL position
    # (force codons for the two downstream amino acids: aa[cp] and aa[cp+1])
    forced_codons_by_gene: Dict[str, Dict[int, str]] = {}
    for gene_id, aa in aa_records_patched:
        aa = aa.strip()
        cuts = cut_points_by_gene[gene_id]
        forced: Dict[int, str] = {}

        for pos_index, cp in enumerate(cuts, start=1):
            oh = assign[pos_index]
            choices = per_gene_choices[gene_id][pos_index]
            if oh not in choices:
                raise RuntimeError(
                    f"Internal error: chosen overhang {oh} not feasible for gene {gene_id} pos {pos_index}"
                )
            choice = choices[oh]
            # Force codons for the two downstream amino acids
            forced[int(cp)] = choice.left_codon
            forced[int(cp + 1)] = choice.right_codon

        forced_codons_by_gene[gene_id] = forced

    return assign, forced_codons_by_gene, aa_records_patched


def run_pipeline(cfg: dict, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)

    seed_codonopt = effective_seed(cfg["codonopt"].get("seed"), cfg["defaults"]["default_seed"])
    seed_reasm = effective_seed(cfg["reassembly"].get("seed"), cfg["defaults"]["default_seed"])

    aa_records_raw = _load_aa_inputs(cfg)

    # Precompute cut points per gene
    cut_points_by_gene: Dict[str, List[int]] = {}
    for gene_id, aa in aa_records_raw:
        cut_points_by_gene[gene_id] = _get_cut_points_aa(cfg, gene_id)

    # Build fragments_aa from AA records (NOTE: these may be patched by consensus during overhang selection)
    # We'll regenerate fragments_aa after overhang selection to ensure we use aa_records_patched consistently.
    # First, compute n_fragments from any gene's cuts:
    any_gene = aa_records_raw[0][0]
    n_fragments = len(cut_points_by_gene[any_gene]) + 1

    # Choose GLOBAL overhangs per junction position (0..k) AND forced codons per gene/AA-index
    overhang_by_pos, forced_codons_by_gene, aa_records = _select_global_overhangs_and_forced_codons(
        cfg=cfg,
        aa_records=aa_records_raw,
        n_fragments=n_fragments,
        cut_points_by_gene=cut_points_by_gene,
    )

    # Now build AA fragments using the possibly patched AA records
    fragments_aa: List[Dict] = []
    for gene_id, aa in aa_records:
        cut_points = cut_points_by_gene[gene_id]
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

    # ---- Overhang audit output + crosstalk threshold ----
    from frogger.overhangs import crosstalk_rows

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
                "fragment_id": frag_id,
                "gene_id": gene_id,
                "frag_index": f["frag_index"],
                "start": f["aa_start"],
                "end": f["aa_end"],
                "core_seq": patched_cds,
                "aa_seq": f["aa_seq"],
                "forced_codons_local": forced_local,
                "aa_len": len(f["aa_seq"]),
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
            return aa

        def _codonopt_runner(fragment_dict: dict):
            new_nt, _, _ = run_codonopt_for_fragment(
                codonopt_cfg=c_cfg,
                outdir=outdir,
                seed=int(seed_codonopt),
                fragment_id=str(fragment_dict["fragment_id"]) + "__repair",
                aa_seq=str(fragment_dict["aa_seq"]),
            )
            return new_nt, ""

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
            repair_path = outdir / cfg.get("outputs", {}).get("repair_report_tsv", "repair_report.tsv")
            write_tsv(repair_path, repair_rows)
    # --------------------------------------------------------------------------

    # Primers / barcodes
    from frogger.primers import build_avoid_kmers, generate_primers_by_pos, primers_from_barcodes

    primers_cfg = cfg.get("primers", {}) or {}
    mode = (primers_cfg.get("mode") or "barcodes").lower()

    if mode == "fixed":
        primers_by_pos = primers_cfg.get("fixed_by_pos") or {}

    elif mode == "barcodes":
        # (6) Default to barcode 901 (forward) and 902 (reverse-complement)
        b = primers_cfg.get("barcodes", {}) or {}
        barcode_fasta = b.get("fasta") or b.get("fasta_path") or primers_cfg.get("barcodes_fasta")
        if not barcode_fasta:
            raise ValueError(
                "primers.mode=barcodes requires primers.barcodes.fasta (path to FASTA with barcode IDs)."
            )
        primers_by_pos = primers_from_barcodes(
            n_positions=n_fragments,
            barcode_fasta=Path(barcode_fasta),
            forward_id=str(b.get("forward_id", "901")),
            reverse_id=str(b.get("reverse_id", "902")),
            forward_seq=b.get("forward_seq"),
            reverse_seq=b.get("reverse_seq"),
            reverse_is_revcomp=bool(b.get("reverse_is_revcomp", True)),
        )

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

    # Build cloning sequences (primer + TypeIIS + OH + core + revcomp(OH) + TypeIIS + primer)
    # (3) The revcomp of the GG overhang is appended to the 3' end of each block in gg_adapters.build_cloning_seq.
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
