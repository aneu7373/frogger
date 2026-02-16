from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

from frogger.common.codon_table import load_codon_usage_xlsx
from frogger.common.gg.splitters import split_aa_by_cut_points
from frogger.gg_adapters import build_cloning_sequences_for_fragments
from frogger.junction_overhangs import enumerate_boundary_overhangs
from frogger.overhangs import load_pairing_matrix_xlsx, choose_overhangs_by_position, OverhangFilters
from frogger.primers import build_avoid_kmers, generate_primers_by_pos


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


def select_global_overhangs_and_forced_codons(
    cfg: dict,
    aa_records: List[Tuple[str, str]],
    n_fragments: int,
    cut_points_by_gene: Dict[str, List[int]],
) -> Tuple[Dict[int, str], Dict[str, Dict[int, str]]]:
    """
    Exactly the shuffle logic, but usable by both shuffle and mutate.

    Returns:
      - overhang_by_pos: {pos_index(0..k): 'ACGT'}
      - forced_codons_by_gene: {gene_id: {aa_index: 'ATG', ...}}
    """
    gg = cfg.get("golden_gate", {}) or {}
    matrix_path = gg.get("overhang_matrix_xlsx")
    if not matrix_path:
        raise ValueError("golden_gate.overhang_matrix_xlsx is required.")
    df = load_pairing_matrix_xlsx(Path(matrix_path), sheet_name=gg.get("overhang_matrix_sheet"))

    c_cfg = cfg.get("codonopt", {}) or {}
    codon_table_path = c_cfg.get("codon_table")
    if not codon_table_path:
        raise ValueError("codonopt.codon_table is required (used to enumerate synonymous codons near junctions).")
    codon_map = load_codon_usage_xlsx(Path(codon_table_path), sheet_name=c_cfg.get("codon_table_sheet"))
    avoid_codons = c_cfg.get("avoid_codons", []) or []

    per_gene_choices: Dict[str, Dict[int, Dict[str, object]]] = {}

    for gene_id, aa in aa_records:
        aa = aa.strip()
        cuts = cut_points_by_gene[gene_id]
        if len(aa) < 2:
            raise ValueError(f"Sequence '{gene_id}' too short for terminus overhang search (<2 aa).")

        boundaries: List[Tuple[int, str, str]] = []
        boundaries.append((0, aa[0], aa[1]))
        for idx, cp in enumerate(cuts, start=1):
            if cp <= 0 or cp >= len(aa):
                raise ValueError(f"Cut point {cp} out of bounds for '{gene_id}' length {len(aa)}")
            boundaries.append((idx, aa[cp - 1], aa[cp]))
        boundaries.append((len(cuts) + 1, aa[-2], aa[-1]))

        per_gene_choices[gene_id] = {}
        for pos_index, aL, aR in boundaries:
            cand_map = enumerate_boundary_overhangs(aL, aR, codon_map, avoid_codons=avoid_codons)
            per_gene_choices[gene_id][pos_index] = cand_map

    # Candidate sets per junction position (intersection across genes; for 1 gene this is just its set)
    candidates_by_pos: Dict[int, List[str]] = {}
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

    forced_codons_by_gene: Dict[str, Dict[int, str]] = {}
    for gene_id, aa in aa_records:
        aa = aa.strip()
        cuts = cut_points_by_gene[gene_id]
        forced: Dict[int, str] = {}

        def aa_indices_for_pos(pos_index: int) -> Tuple[int, int]:
            if pos_index == 0:
                return 0, 1
            if pos_index == n_fragments:
                return len(aa) - 2, len(aa) - 1
            cp = cuts[pos_index - 1]
            return cp - 1, cp

        for pos_index in range(0, n_fragments + 1):
            oh = assign[pos_index]
            choices = per_gene_choices[gene_id][pos_index]
            choice = choices[oh]
            iL, iR = aa_indices_for_pos(pos_index)
            forced[iL] = choice.left_codon
            forced[iR] = choice.right_codon

        forced_codons_by_gene[gene_id] = forced

    return assign, forced_codons_by_gene


def build_primers_by_pos_for_fragments(cfg: dict, fragments_cds: List[Dict], n_fragments: int) -> Dict[int, Dict[str, str]]:
    primers_cfg = cfg.get("primers", {}) or {}
    mode = (primers_cfg.get("mode") or "fixed").lower()

    if mode == "fixed":
        return primers_cfg.get("fixed_by_pos") or {}

    if mode == "generate":
        g = primers_cfg.get("generate", {}) or {}
        k = int(g.get("max_internal_match", 12))
        avoid = build_avoid_kmers([f["core_seq"] for f in fragments_cds], k=k)
        return generate_primers_by_pos(
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

    raise ValueError(f"Unknown primers.mode: {mode}")


def compute_overhangs_by_frag_pos(overhang_by_pos: Dict[int, str], n_fragments: int) -> Dict[int, Dict[str, str]]:
    out: Dict[int, Dict[str, str]] = {}
    for i in range(1, n_fragments + 1):
        out[i] = {"left": overhang_by_pos[i - 1], "right": overhang_by_pos[i]}
    return out
