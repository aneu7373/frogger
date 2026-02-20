from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

from frogger.common.codon_table import load_codon_usage_xlsx
from frogger.gg_adapters import build_cloning_sequences_for_fragments
from frogger.junction_overhangs import enumerate_boundary_overhangs
from frogger.overhangs import load_pairing_matrix_xlsx, choose_overhangs_by_position, OverhangFilters
from frogger.primers import build_avoid_kmers, generate_primers_by_pos, primers_from_barcodes


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


def _consensus_downstream_aas(
    aa_records: List[Tuple[str, str]],
    cut_points_by_gene: Dict[str, List[int]],
) -> Dict[int, Tuple[str, str]]:
    """
    For each internal junction (pos_index 1..len(cuts)), choose the most common amino acids
    for the FIRST TWO amino acids AFTER the cut site.

    For a cut point cp (between cp-1 | cp), downstream amino acids are:
      aa1 = seq[cp]
      aa2 = seq[cp+1]
    """
    any_gene = aa_records[0][0]
    n_internal = len(cut_points_by_gene[any_gene])

    consensus: Dict[int, Tuple[str, str]] = {}
    for pos_index in range(1, n_internal + 1):
        a1 = Counter()
        a2 = Counter()
        for gene_id, aa in aa_records:
            aa = aa.strip()
            cuts = cut_points_by_gene[gene_id]
            cp = int(cuts[pos_index - 1])
            if cp <= 0 or cp + 1 >= len(aa):
                raise ValueError(
                    f"Cut point {cp} invalid for downstream-AA window in '{gene_id}' length {len(aa)}. "
                    "Cut points must allow two amino acids after the cut (cp <= len(aa)-2)."
                )
            a1[aa[cp]] += 1
            a2[aa[cp + 1]] += 1

        aa1 = sorted(a1.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        aa2 = sorted(a2.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        consensus[pos_index] = (aa1, aa2)
    return consensus


def apply_consensus_aas_to_records(
    aa_records: List[Tuple[str, str]],
    cut_points_by_gene: Dict[str, List[int]],
) -> List[Tuple[str, str]]:
    """
    Apply consensus downstream AAs (first two AAs after each cut) to each AA sequence.

    Requested behavior for shuffle:
      - If amino acids are different at a cut site between genes, use the most common in the set.
      - It is OK if this changes the gene sequence (proteins will be shuffled anyway).

    Returns a NEW aa_records list with sequences patched at [cp] and [cp+1] for each cut.
    """
    if len(aa_records) <= 1:
        return aa_records

    consensus = _consensus_downstream_aas(aa_records, cut_points_by_gene)

    patched: List[Tuple[str, str]] = []
    for gene_id, aa in aa_records:
        aa = aa.strip()
        arr = list(aa)
        cuts = cut_points_by_gene[gene_id]
        for pos_index, (aa1, aa2) in consensus.items():
            cp = int(cuts[pos_index - 1])
            arr[cp] = aa1
            arr[cp + 1] = aa2
        patched.append((gene_id, "".join(arr)))
    return patched


def select_global_overhangs_and_forced_codons(
    cfg: dict,
    aa_records: List[Tuple[str, str]],
    n_fragments: int,
    cut_points_by_gene: Dict[str, List[int]],
) -> Tuple[Dict[int, str], Dict[str, Dict[int, str]], List[Tuple[str, str]]]:
    """
    Unified overhang + forced-codon selection for shuffle and mutate.

    Updates implemented:
      1) Overhang search window = FIRST TWO amino acids AFTER cut.
      2) 5' and 3' terminal overhangs are REQUIRED and not chosen by matrix search.
      3) Matrix compatibility scoring includes reverse complements (handled in overhangs.py).
      4) Codons used during overhang enumeration are restricted to top-K most frequent codons (default K=2).
      5) In multi-gene shuffle, if AAs differ at a cut site, use the most common (and patch AA sequences).

    Returns:
      - overhang_by_pos: {pos_index 0..n_fragments}
      - forced_codons_by_gene: {gene_id: {aa_index: codon}} (internal junctions only; indices are global AA indices)
      - aa_records_patched: patched AA records (same as input for single gene)
    """
    if not aa_records:
        raise ValueError("No AA records provided.")

    term5, term3 = _require_terminal_overhangs(cfg)

    gg = cfg.get("golden_gate", {}) or {}
    matrix_path = gg.get("overhang_matrix_xlsx")
    if not matrix_path:
        raise ValueError("golden_gate.overhang_matrix_xlsx is required.")
    df = load_pairing_matrix_xlsx(Path(matrix_path), sheet_name=gg.get("overhang_matrix_sheet"))

    aa_records_patched = apply_consensus_aas_to_records(aa_records, cut_points_by_gene)

    c_cfg = cfg.get("codonopt", {}) or {}
    codon_table_path = c_cfg.get("codon_table")
    if not codon_table_path:
        raise ValueError("codonopt.codon_table is required (used to enumerate synonymous codons near junctions).")
    codon_map = load_codon_usage_xlsx(Path(codon_table_path), sheet_name=c_cfg.get("codon_table_sheet"))
    avoid_codons = c_cfg.get("avoid_codons", []) or []

    search_cfg = gg.get("overhang_search", {}) or {}
    top_k = int(search_cfg.get("top_codons_per_aa", 2))

    per_gene_choices: Dict[str, Dict[int, Dict[str, object]]] = {}
    for gene_id, aa in aa_records_patched:
        aa = aa.strip()
        cuts = cut_points_by_gene[gene_id]
        if len(cuts) != n_fragments - 1:
            raise ValueError(
                f"Gene '{gene_id}' has {len(cuts)} cut points but n_fragments={n_fragments} "
                f"(expected {n_fragments-1} cut points)."
            )

        for cp in cuts:
            if cp <= 0 or cp + 1 >= len(aa):
                raise ValueError(
                    f"Cut point {cp} invalid for downstream-AA window in '{gene_id}' length {len(aa)}. "
                    "Cut points must allow two amino acids after the cut (cp <= len(aa)-2)."
                )

        per_gene_choices[gene_id] = {}
        for pos_index, cp in enumerate(cuts, start=1):
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

    candidates_by_pos: Dict[int, List[str]] = {}
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

    flt_cfg = (gg.get("overhang_filters", {}) or {})
    flt = OverhangFilters(
        disallow_palindromes=bool(flt_cfg.get("disallow_palindromes", True)),
        disallow_revcomp_self=bool(flt_cfg.get("disallow_revcomp_self", True)),
        min_gc=float(flt_cfg.get("min_gc", 0.25)),
        max_gc=float(flt_cfg.get("max_gc", 0.75)),
        max_homopolymer=int(flt_cfg.get("max_homopolymer", 3)),
        # Optional: explicitly ban certain 4-mers (helps eliminate matrix “landmines” like CATT)
        disallow_overhangs=tuple(str(x).upper().strip() for x in (flt_cfg.get("disallow_overhangs", []) or [])),
    )

    # Include fixed terminal overhangs inside the optimizer objective so internal picks cannot
    # massively cross-talk with AATG/GCTT (or collide via reverse complements).
    # Optionally ignore terminal↔terminal cross-talk in the worst-score objective, since those
    # ends generally ligate to a vector backbone rather than to each other.
    ignore_tt = bool(gg.get("ignore_terminal_terminal_crosstalk", True))

    internal_assign, worst = choose_overhangs_by_position(
        df=df,
        candidates_by_pos=candidates_by_pos,
        filters=flt,
        beam_width=int(gg.get("beam_width", 200)),
        fixed_overhangs=[term5, term3],
        include_fixed_fixed_in_worst=(not ignore_tt),
    )

    assign: Dict[int, str] = {0: term5, n_fragments: term3}
    for pos_index, oh in internal_assign.items():
        assign[int(pos_index)] = str(oh).upper()

    gg["_selected_overhangs_by_pos"] = assign
    gg["_selected_overhangs_worst_score"] = worst
    cfg["golden_gate"] = gg

    forced_codons_by_gene: Dict[str, Dict[int, str]] = {}
    for gene_id, aa in aa_records_patched:
        aa = aa.strip()
        cuts = cut_points_by_gene[gene_id]
        forced: Dict[int, str] = {}
        for pos_index, cp in enumerate(cuts, start=1):
            oh = assign[pos_index]
            choices = per_gene_choices[gene_id][pos_index]
            choice = choices[oh]
            forced[int(cp)] = choice.left_codon
            forced[int(cp + 1)] = choice.right_codon
        forced_codons_by_gene[gene_id] = forced

    return assign, forced_codons_by_gene, aa_records_patched


def build_primers_by_pos_for_fragments(cfg: dict, fragments_cds: List[Dict], n_fragments: int) -> Dict[int, Dict[str, str]]:
    primers_cfg = cfg.get("primers", {}) or {}
    mode = (primers_cfg.get("mode") or "barcodes").lower()

    if mode == "fixed":
        return primers_cfg.get("fixed_by_pos") or {}

    if mode == "barcodes":
        b = primers_cfg.get("barcodes", {}) or {}
        barcode_fasta = b.get("fasta") or b.get("fasta_path") or primers_cfg.get("barcodes_fasta")
        if not barcode_fasta:
            raise ValueError(
                "primers.mode=barcodes requires primers.barcodes.fasta (path to FASTA with barcode IDs)."
            )
        forward_id = str(b.get("forward_id", "901"))
        reverse_id = str(b.get("reverse_id", "902"))
        forward_seq = b.get("forward_seq")
        reverse_seq = b.get("reverse_seq")
        reverse_is_revcomp = bool(b.get("reverse_is_revcomp", True))
        return primers_from_barcodes(
            n_positions=n_fragments,
            barcode_fasta=Path(barcode_fasta),
            forward_id=forward_id,
            reverse_id=reverse_id,
            forward_seq=forward_seq,
            reverse_seq=reverse_seq,
            reverse_is_revcomp=reverse_is_revcomp,
        )

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
    """
    Convert junction overhang assignment (pos 0..n_fragments) into per-fragment left/right overhangs.

    NOTE:
      - build_cloning_seq() appends revcomp(right_overhang) to the 3' end of each block.
      - Therefore we store right_overhang in forward orientation here (the designed overhang label).
    """
    out: Dict[int, Dict[str, str]] = {}
    for i in range(1, n_fragments + 1):
        out[i] = {"left": overhang_by_pos[i - 1], "right": overhang_by_pos[i]}
    return out
