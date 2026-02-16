from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd


DNA = str


def revcomp(seq: DNA) -> DNA:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def gc_fraction(seq: DNA) -> float:
    s = seq.upper()
    if not s:
        return 0.0
    return (s.count("G") + s.count("C")) / len(s)


@dataclass(frozen=True)
class OverhangFilters:
    disallow_palindromes: bool = True
    disallow_revcomp_self: bool = True
    min_gc: float = 0.25
    max_gc: float = 0.75
    max_homopolymer: int = 3  # blocks AAAA/TTTT/CCCC/GGGG by default


def _max_homopolymer_run(seq: DNA) -> int:
    best = 1
    cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


def load_pairing_matrix_xlsx(path: Path, sheet_name: Optional[str] = None) -> pd.DataFrame:
    """
    Loads a 256x256 4-mer pairing matrix from an Excel file.
    Expected format: index column = overhangs, columns = overhangs, values = numeric score (lower is better).
    """
    xls = pd.ExcelFile(path)
    sheet = sheet_name or xls.sheet_names[0]
    df = pd.read_excel(path, sheet_name=sheet, index_col=0)
    df.index = df.index.astype(str).str.upper()
    df.columns = df.columns.astype(str).str.upper()
    return df


def filter_overhangs(candidates: Iterable[DNA], flt: OverhangFilters) -> List[DNA]:
    out: List[DNA] = []
    for oh in candidates:
        oh = str(oh).upper().strip()
        if len(oh) != 4 or any(b not in "ACGT" for b in oh):
            continue
        if flt.disallow_palindromes and oh == oh[::-1]:
            continue
        if flt.disallow_revcomp_self and oh == revcomp(oh):
            continue
        gc = gc_fraction(oh)
        if gc < flt.min_gc or gc > flt.max_gc:
            continue
        if _max_homopolymer_run(oh) > flt.max_homopolymer:
            continue
        out.append(oh)
    return sorted(set(out))


def symmetric_pair_score(df: pd.DataFrame, a: DNA, b: DNA) -> float:
    """
    Matrix may be non-symmetric; use the worse of the two directions.
    """
    a = a.upper()
    b = b.upper()
    v1 = float(df.loc[a, b])
    v2 = float(df.loc[b, a])
    return max(v1, v2)

def crosstalk_rows(df: pd.DataFrame, assign: Dict[int, DNA]) -> List[Dict[str, object]]:
    """
    Build pairwise cross-talk rows for a chosen overhang assignment.
    Uses the same symmetric scoring as the search (max of both directions).
    """
    items = sorted(assign.items(), key=lambda kv: kv[0])
    overhangs = [oh for _pos, oh in items]
    rows: List[Dict[str, object]] = []
    for i in range(len(overhangs)):
        for j in range(i + 1, len(overhangs)):
            a = overhangs[i]
            b = overhangs[j]
            rows.append(
                {
                    "overhang_a": a,
                    "overhang_b": b,
                    "score": float(symmetric_pair_score(df, a, b)),
                }
            )
    rows.sort(key=lambda r: r["score"], reverse=True)
    return rows

def choose_overhangs_by_position(
    df: pd.DataFrame,
    candidates_by_pos: Dict[int, Sequence[DNA]],
    filters: Optional[OverhangFilters] = None,
    beam_width: int = 200,
) -> Tuple[Dict[int, DNA], float]:
    """
    Choose one overhang per position from candidates_by_pos, minimizing the maximum symmetric cross-talk score
    among all chosen overhangs. Enforces all chosen overhangs are distinct.

    candidates_by_pos: {pos_index: [4mers...]} pos_index is in assembly order (0..k).
    Returns ({pos: overhang}, worst_pair_score)
    """
    flt = filters or OverhangFilters()
    # Pre-filter candidates per position
    filtered: Dict[int, List[DNA]] = {}
    for pos, cands in candidates_by_pos.items():
        filtered[pos] = filter_overhangs(cands, flt)
        if not filtered[pos]:
            raise ValueError(f"No feasible overhang candidates for position {pos} after filtering.")

    positions = sorted(filtered.keys())

    # Precompute pair score cache for all overhangs that appear
    all_oh = sorted(set([oh for pos in positions for oh in filtered[pos]]))
    pair_cache: Dict[Tuple[DNA, DNA], float] = {}
    for i, a in enumerate(all_oh):
        for b in all_oh[i + 1 :]:
            pair_cache[(a, b)] = symmetric_pair_score(df, a, b)

    def pair(a: DNA, b: DNA) -> float:
        if a == b:
            return float("inf")
        x, y = (a, b) if a < b else (b, a)
        return pair_cache[(x, y)]

    # Beam states are (assignment_dict, chosen_list, worst_so_far)
    beam = [({}, [], 0.0)]
    for pos in positions:
        new_beam = []
        for assign, chosen, worst in beam:
            chosen_set = set(chosen)
            for cand in filtered[pos]:
                if cand in chosen_set:
                    continue
                new_worst = worst
                for prev in chosen:
                    new_worst = max(new_worst, pair(prev, cand))
                new_assign = dict(assign)
                new_assign[pos] = cand
                new_beam.append((new_assign, chosen + [cand], new_worst))
        # Keep best states by worst score
        new_beam.sort(key=lambda x: x[2])
        beam = new_beam[: max(1, int(beam_width))]

        if not beam:
            raise RuntimeError(f"Overhang search failed at position {pos} (beam exhausted).")

    best_assign, _, best_worst = beam[0]
    return best_assign, float(best_worst)
