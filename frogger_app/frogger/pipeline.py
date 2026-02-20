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

def _translate_nt_strict(nt: str, *, context: str) -> str:
    nt = str(nt)
    if len(nt) % 3 != 0:
        raise ValueError(f"Refusing to translate non-multiple-of-3 nt (len={len(nt)}) in {context}")
    return str(Seq(nt).translate(to_stop=False))

def _split_cds_fix_a(
    cds_full: str,
    *,
    aa_cut_points: List[int],
    window_offset_by_pos: Dict[int, int],
) -> List[Dict[str, int | str]]:
    """
    FIX A splitter.

    We split a *full-length* CDS into Golden-Gate fragment cores such that each internal junction
    is a 4-nt overlap derived from the native CDS sequence (not appended).
    - For junction pos_index=1..k-1 at AA cut point cp (0-based AA index),
      compute the junction's nucleotide index:
          j = 3*cp + window_offset_by_pos[pos_index]
      Then:
        - left fragment includes cds_full[... : j+4]
        - right fragment starts at cds_full[j : ...]
      so the 4-nt region cds_full[j:j+4] is shared (the physical sticky region).
    """
    seq = (cds_full or "").strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
    if len(seq) % 3 != 0:
        raise ValueError(f"Full CDS must be multiple-of-3 (len={len(seq)}).")

    cps_aa = [int(x) for x in (aa_cut_points or [])]
    if sorted(cps_aa) != cps_aa:
        raise ValueError("aa_cut_points must be sorted ascending.")
    if len(set(cps_aa)) != len(cps_aa):
        raise ValueError("aa_cut_points must not contain duplicates.")

    aa_len = len(seq) // 3
    for cp in cps_aa:
        if cp <= 0 or cp >= aa_len:
            raise ValueError(f"Invalid AA cut point {cp} for aa_len={aa_len}. Must be within (0, aa_len).")

    # Internal junction nt indices
    j_nt: List[int] = []
    for pos_index, cp in enumerate(cps_aa, start=1):
        off = int(window_offset_by_pos.get(pos_index, 0))
        if off not in (0, 1, 2):
            raise ValueError(f"window_offset_by_pos[{pos_index}] must be 0/1/2, got {off}")
        j = 3 * int(cp) + off
        if j < 0 or j + 4 > len(seq):
            raise ValueError(f"Computed junction nt index out of range: pos_index={pos_index} j={j} len={len(seq)}")
        j_nt.append(j)

    # Fragment AA boundaries (frame-aligned)
    aa_points = [0] + cps_aa + [aa_len]

    out: List[Dict[str, int | str]] = []
    n_frags = len(aa_points) - 1
    for i in range(1, n_frags + 1):
        aa_a = aa_points[i - 1]
        aa_b = aa_points[i]

        nt_start = 0 if i == 1 else j_nt[i - 2]
        nt_end = len(seq) if i == n_frags else (j_nt[i - 1] + 4)

        core = seq[nt_start:nt_end]
        out.append(
            {
                "frag_index": i,
                "aa_start0": aa_a,
                "aa_end0": aa_b,
                "nt_start0": nt_start,
                "nt_end0": nt_end,
                "core_seq": core,
            }
        )
    return out

def _reassemble_fix_a_cores(frs_sorted: List[Dict]) -> str:
    """Reassemble Fix-A cores back into a full CDS (drop the 4-nt overlap from all but the first)."""
    if not frs_sorted:
        return ""
    return str(frs_sorted[0]["core_seq"]) + "".join(str(f["core_seq"])[4:] for f in frs_sorted[1:])


def _repair_full_gene_cds(
    *,
    gene_id: str,
    aa: str,
    cds_full: str,
    c_cfg: dict,
    outdir: Path,
    seed_codonopt: int,
    repair_cfg: dict,
) -> Tuple[str, List[dict]]:
    """
    Run existing repair loop on a SINGLE full-gene CDS (multiple-of-3), not on overlapped fragments.
    Returns (cds_repaired, repair_rows).
    """
    from frogger.common.repair import repair_fragments_loop

    qc_cfg = (repair_cfg.get("qc", {}) or {})
    max_rounds = int(repair_cfg.get("max_rounds", 1))

    # Use repair_fragments_loop on a single "fragment" representing the full gene.
    frag = {
        "fragment_id": f"{gene_id}__FULL",
        "gene_id": gene_id,
        "frag_index": 1,
        "start": 0,
        "end": len(aa),
        "core_seq": cds_full,
        "aa_seq": aa,
        "forced_codons_local": {},  # not used for full-gene repair
        "aa_len": len(aa),
    }

    def _translate_nt(nt: str) -> str:
        # Be STRICT: if repair ever creates non-multiple-of-3 intermediates, fail immediately
        # rather than letting Biopython emit warnings.
        return _translate_nt_strict(nt, context=f"{gene_id}:repair_full_gene")


    def _codonopt_runner(fragment_dict: dict):
        new_nt, _, _ = run_codonopt_for_fragment(
            codonopt_cfg=c_cfg,
            outdir=outdir,
            seed=int(seed_codonopt),
            fragment_id=str(fragment_dict["fragment_id"]) + "__repair",
            aa_seq=str(fragment_dict["aa_seq"]),
        )
        return new_nt, ""

    repaired_list, repair_rows = repair_fragments_loop(
        [frag],
        qc_cfg=qc_cfg,
        enable=True,
        max_rounds=max_rounds,
        codonopt_fragment_runner=_codonopt_runner,
        aa_translate=_translate_nt,
        patch_forced_codons=_patch_forced_codons,
    )

    if not repaired_list:
        raise RuntimeError("repair_fragments_loop returned empty list for full gene.")

    cds_repaired = str(repaired_list[0]["core_seq"])
    return cds_repaired, (repair_rows or [])

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

def _group_frags_by_gene_id(frags: List[dict]) -> Dict[str, List[dict]]:
    groups: Dict[str, List[dict]] = {}
    for f in frags:
        gid = str(f.get("gene_id") or f.get("parent_id") or f.get("id") or "UNKNOWN")
        groups.setdefault(gid, []).append(f)
    return groups

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
        fixed_overhangs=[term5, term3],
        # Don't count term5↔term3 as part of the minimax objective (they never meet in a linear assembly)
        include_fixed_fixed_in_worst=False,
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
    #
    # FIX A support:
    #   - also return window_offset_by_pos so we can split CDS at nucleotide indices that embed the overhang
    #   - enforce that all genes use the SAME window_offset for a given junction position (critical for shuffle pools)

    window_offset_by_pos: Dict[int, int] = {}
    forced_codons_by_gene: Dict[str, Dict[int, str]] = {}

    # Pick a canonical reference gene for offsets/codons (deterministic)
    ref_gene_id = aa_records_patched[0][0]

    # Precompute canonical (pos_index -> JunctionChoice) from ref gene
    canonical_choice_by_pos: Dict[int, object] = {}
    ref_cuts = cut_points_by_gene[ref_gene_id]
    for pos_index, cp in enumerate(ref_cuts, start=1):
        oh = assign[pos_index]
        ref_choices = per_gene_choices[ref_gene_id][pos_index]
        if oh not in ref_choices:
            raise RuntimeError(
                f"Internal error: chosen overhang {oh} not feasible for reference gene {ref_gene_id} pos {pos_index}"
            )
        canonical_choice_by_pos[pos_index] = ref_choices[oh]
        window_offset_by_pos[pos_index] = int(ref_choices[oh].window_offset)

    # Apply canonical codons + validate offsets across all genes
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

            ch = choices[oh]
            canon = canonical_choice_by_pos[pos_index]

            # Enforce identical window_offset across genes (otherwise pooled fragments can create mixed codons)
            if int(ch.window_offset) != int(canon.window_offset):
                raise ValueError(
                    "Fix A requires consistent window_offset across all genes for a given junction.\n"
                    f"  pos_index={pos_index}\n"
                    f"  chosen_overhang={oh}\n"
                    f"  ref_gene={ref_gene_id} window_offset={canon.window_offset}\n"
                    f"  this_gene={gene_id} window_offset={ch.window_offset}\n"
                    "Try increasing golden_gate.overhang_search.top_codons_per_aa, or change cut points."
                )

            # Force CANONICAL codons for the two downstream amino acids
            forced[int(cp)] = canon.left_codon
            forced[int(cp + 1)] = canon.right_codon

        forced_codons_by_gene[gene_id] = forced

    # stash for audit/debug
    gg["_selected_window_offset_by_pos"] = window_offset_by_pos
    cfg["golden_gate"] = gg

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

    gg = cfg.get("golden_gate", {}) or {}
    window_offset_by_pos = gg.get("_selected_window_offset_by_pos")
    if not window_offset_by_pos:
        raise RuntimeError("Missing golden_gate._selected_window_offset_by_pos after overhang selection.")

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

    from frogger.overhangs import rc_aware_pair_matrix, desired_pair_score, revcomp

    selected_overhangs_in_order = [oh for _pos, oh in sorted(assign.items(), key=lambda kv: kv[0])]

    # (A) Full matrix including reverse complements as explicit rows/cols
    mat_rc = rc_aware_pair_matrix(df_matrix, selected_overhangs_in_order, include_reverse_complements=True)
    (outdir / "overhang_crosstalk_matrix_with_rc.csv").write_text(mat_rc.to_csv(index=True))

    # (B) Intended pair strengths: a vs rc(a)
    rows = []
    for oh in selected_overhangs_in_order:
        rows.append(
            {
                "overhang": oh,
                "revcomp": revcomp(oh),
                "desired_pair_score": float(desired_pair_score(df_matrix, oh)),
            }
        )
    write_tsv(outdir / "overhang_desired_pairs.tsv", rows)
    # Map junction overhangs into fragment-end overhangs
    overhangs_by_frag_pos: Dict[int, Dict[str, str]] = {}
    for i in range(1, n_fragments + 1):
        overhangs_by_frag_pos[i] = {"left": assign[i - 1], "right": assign[i]}

    # Codonopt per fragment, then patch forced codons for AA positions that fall inside this fragment
    # Codonopt FULL gene (per gene), patch forced codons globally, then split with FIX A overlap logic.
    fragments_cds: List[Dict] = []
    c_cfg = cfg["codonopt"]

    # Build a quick lookup of AA fragment boundaries per gene/frag_index (for headers + reporting)
    aa_bounds_by_gene_and_pos: Dict[tuple[str, int], tuple[int, int, str]] = {}
    for f in fragments_aa:
        aa_bounds_by_gene_and_pos[(str(f["gene_id"]), int(f["frag_index"]))] = (
            int(f["aa_start"]),
            int(f["aa_end"]),
            str(f["aa_seq"]),
        )

    for gene_id, aa in aa_records:
        gene_id = str(gene_id)
        aa = aa.strip()

        # 1) codonopt the whole ORF
        gene_frag_id = f"{gene_id}__FULL"
        cds_full, _, _ = run_codonopt_for_fragment(
            codonopt_cfg=c_cfg,
            outdir=outdir,
            seed=int(seed_codonopt),
            fragment_id=gene_frag_id,
            aa_seq=aa,
        )

        # 2) patch forced codons at junction-downstream AA positions (global AA indices)
        forced_global = forced_codons_by_gene.get(gene_id, {})
        cds_patched = _patch_forced_codons(cds_full, aa_len=len(aa), forced_by_local_aa=forced_global)

        # 3) sanity: translation of full patched CDS must match AA
        if _translate_nt_strict(cds_patched, context=f"{gene_id}:cds_patched") != aa:
            raise RuntimeError(
                f"Forced-codon patching changed translation for {gene_id} FULL gene. "
                f"This indicates an internal bug in codon forcing."
            )

        # 4) FIX A split: overlapping 4-nt junctions embedded in the native CDS
        aa_cut_points = cut_points_by_gene[gene_id]
        frag_dicts = _split_cds_fix_a(
            cds_patched,
            aa_cut_points=aa_cut_points,
            window_offset_by_pos=window_offset_by_pos,
        )

        # 5) build fragment records in the shape expected downstream
        for d in frag_dicts:
            pos = int(d["frag_index"])
            aa_start, aa_end, frag_aa = aa_bounds_by_gene_and_pos[(gene_id, pos)]
            frag_id = f"{gene_id}__pos{pos:02d}__aa{aa_start}-{aa_end}"

            fragments_cds.append(
                {
                    "fragment_id": frag_id,
                    "gene_id": gene_id,
                    "frag_index": pos,
                    "start": aa_start,
                    "end": aa_end,
                    "core_seq": str(d["core_seq"]),
                    "aa_seq": frag_aa,
                    "forced_codons_local": {},  # no longer meaningful per-fragment under Fix A
                    "aa_len": len(frag_aa),
                    "nt_start0": int(d["nt_start0"]),
                    "nt_end0": int(d["nt_end0"]),
                }
            )

        # 6) sanity: reassemble cores by removing the 4-nt overlap from all but the first
        #    (this must exactly reproduce cds_patched)
        frs = [x for x in fragments_cds if x["gene_id"] == gene_id]
        frs.sort(key=lambda z: int(z["frag_index"]))
        if not frs:
            raise RuntimeError(f"No fragments created for gene {gene_id}")
        cds_reasm = frs[0]["core_seq"] + "".join(str(f["core_seq"])[4:] for f in frs[1:])
        if cds_reasm != cds_patched:
            raise RuntimeError(
                f"FIX A split/reassembly mismatch for gene {gene_id}:\n"
                f"  expected_len={len(cds_patched)} got_len={len(cds_reasm)}\n"
                f"  This indicates incorrect window_offset_by_pos usage or splitter logic."
            )

    # --- QC + repair pass (Fix A compatible): run repair on FULL gene CDS, then re-split -------------
    repair_cfg = (cfg.get("repair", {}) or {})
    repair_rows_all: List[dict] = []

    if bool(repair_cfg.get("enable", False)):
        # Build a quick mapping of gene_id -> aa (patched AA)
        aa_by_gene = {str(g): str(a).strip() for g, a in aa_records}

        # Reconstruct per-gene full CDS from fragments_cds (which currently came from cds_patched),
        # repair it safely (full CDS is multiple-of-3), validate translation unchanged, then re-split.
        new_fragments_cds: List[Dict] = []

        for gene_id, aa in aa_by_gene.items():
            # Grab current fragments in order
            frs = [x for x in fragments_cds if x["gene_id"] == gene_id]
            frs.sort(key=lambda z: int(z["frag_index"]))
            if not frs:
                raise RuntimeError(f"No fragments found for gene {gene_id} during repair.")

            cds_curr = _reassemble_fix_a_cores(frs)

            # Sanity: current full CDS should translate to the gene AA
            if len(cds_curr) % 3 != 0:
                raise RuntimeError(
                    f"Fix-A reassembled full CDS for {gene_id} is not multiple-of-3 (len={len(cds_curr)})."
                )
            if _translate_nt_strict(cds_curr, context=f"{gene_id}:cds_curr") != aa:
                raise RuntimeError(f"Fix-A reassembled full CDS translation mismatch for {gene_id} prior to repair.")

            cds_repaired, repair_rows = _repair_full_gene_cds(
                gene_id=gene_id,
                aa=aa,
                cds_full=cds_curr,
                c_cfg=c_cfg,
                outdir=outdir,
                seed_codonopt=int(seed_codonopt),
                repair_cfg=repair_cfg,
            )
            repair_rows_all.extend(repair_rows)

            # Sanity: repair must preserve translation
            if len(cds_repaired) % 3 != 0:
                raise RuntimeError(f"Repair produced non-multiple-of-3 CDS for {gene_id} (len={len(cds_repaired)}).")
            if _translate_nt_strict(cds_repaired, context=f"{gene_id}:cds_repaired") != aa:
                raise RuntimeError(f"Repair changed translation for {gene_id}. Repair must preserve AA sequence.")

            # Re-split repaired full CDS using Fix A offsets
            aa_cut_points = cut_points_by_gene[gene_id]
            frag_dicts = _split_cds_fix_a(
                cds_repaired,
                aa_cut_points=aa_cut_points,
                window_offset_by_pos=window_offset_by_pos,
            )

            # Rebuild fragment dicts (same shape as before)
            for d in frag_dicts:
                pos = int(d["frag_index"])
                aa_start, aa_end, frag_aa = aa_bounds_by_gene_and_pos[(gene_id, pos)]
                frag_id = f"{gene_id}__pos{pos:02d}__aa{aa_start}-{aa_end}"

                new_fragments_cds.append(
                    {
                        "fragment_id": frag_id,
                        "gene_id": gene_id,
                        "frag_index": pos,
                        "start": aa_start,
                        "end": aa_end,
                        "core_seq": str(d["core_seq"]),
                        "aa_seq": frag_aa,
                        "forced_codons_local": {},
                        "aa_len": len(frag_aa),
                        "nt_start0": int(d["nt_start0"]),
                        "nt_end0": int(d["nt_end0"]),
                    }
                )

            # Final sanity: reassemble equals repaired CDS
            frs2 = [x for x in new_fragments_cds if x["gene_id"] == gene_id]
            frs2.sort(key=lambda z: int(z["frag_index"]))
            if _reassemble_fix_a_cores(frs2) != cds_repaired:
                raise RuntimeError(f"Fix-A re-split/reassembly mismatch for {gene_id} after repair.")

        # Swap in repaired fragments
        fragments_cds = new_fragments_cds

        if repair_rows_all:
            repair_path = outdir / cfg.get("outputs", {}).get("repair_report_tsv", "repair_report.tsv")
            write_tsv(repair_path, repair_rows_all)

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
        n_fragments=n_fragments,
    )
    
    # --- (Optional but recommended) Validate physical GG adjacency per gene -----------------
    # This avoids false failures when your FASTA contains multiple genes (e.g., 3x8 fragments)
    # and ensures we only validate within each gene in positional order.
    from frogger.gg_adapters import validate_physical_golden_gate_chain

    groups = _group_frags_by_gene_id(gg_frags)
    for gid, frs in groups.items():
        frs_sorted = sorted(frs, key=lambda x: int(x.get("frag_index", 0)))
        #validate_physical_golden_gate_chain(frs_sorted, enzyme)
    # ---------------------------------------------------------------------------------------

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
                "core_len_nt": len(f["core_seq"]),
                "cloning_len_nt": len(f["cloning_seq"]),
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