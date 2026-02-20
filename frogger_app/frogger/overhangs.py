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
    # Optional hard ban list (exact 4-mers, uppercase)
    disallow_overhangs: Tuple[str, ...] = ()


def _max_homopolymer_run(seq: DNA) -> int:
    if not seq:
        return 0
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
    banned = set(x.upper() for x in (flt.disallow_overhangs or ()))
    out: List[DNA] = []
    for oh in candidates:
        oh = str(oh).upper().strip()
        if len(oh) != 4 or any(b not in "ACGT" for b in oh):
            continue
        if oh in banned:
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


def rc_aware_pair_score(df: pd.DataFrame, a: DNA, b: DNA) -> float:
    """
    RC-aware scoring:
      - score (a, b) AND also against the reverse-complement representations.
    We conservatively take the worst (maximum) score across:
      (a, b), (a, rc(b)), (rc(a), b), (rc(a), rc(b))
    """
    a = a.upper()
    b = b.upper()
    ra = revcomp(a)
    rb = revcomp(b)
    return max(
        symmetric_pair_score(df, a, b),
        symmetric_pair_score(df, a, rb),
        symmetric_pair_score(df, ra, b),
        symmetric_pair_score(df, ra, rb),
    )

def desired_pair_score(df: pd.DataFrame, a: DNA) -> float:
    """
    Score for the *intended* ligation partner: a vs revcomp(a).

    IMPORTANT: This assumes the matrix values represent "compatibility / ligation propensity"
    where *lower is better* for non-pairs. If your matrix is inverted (higher is better),
    flip the optimization sense accordingly.
    """
    a = a.upper()
    ra = revcomp(a)
    return float(symmetric_pair_score(df, a, ra))


def overhang_list_with_rc(overhangs: Sequence[DNA]) -> List[DNA]:
    """Return [a1, a2, ...] plus any missing reverse-complements, without duplicates."""
    out: List[DNA] = []
    seen = set()
    for x in overhangs:
        x = str(x).upper().strip()
        if x not in seen:
            out.append(x)
            seen.add(x)
        rx = revcomp(x)
        if rx not in seen:
            out.append(rx)
            seen.add(rx)
    return out


def rc_aware_pair_matrix(df: pd.DataFrame, overhangs: Sequence[DNA], *, include_reverse_complements: bool = False) -> pd.DataFrame:
    """
    Square matrix of RC-aware cross-talk scores for the provided overhang list.
    If include_reverse_complements=True, rows/cols include each overhang and its RC explicitly.
    """
    base = [str(x).upper().strip() for x in overhangs]
    ohs = overhang_list_with_rc(base) if include_reverse_complements else base

    mat = pd.DataFrame(index=ohs, columns=ohs, dtype=float)
    for a in ohs:
        for b in ohs:
            if a == b:
                mat.loc[a, b] = 0.0
            else:
                mat.loc[a, b] = float(rc_aware_pair_score(df, a, b))
    return mat


def crosstalk_rows(df: pd.DataFrame, assign: Dict[int, DNA]) -> List[Dict[str, object]]:
    """
    Build pairwise cross-talk rows for a chosen overhang assignment.
    Uses the same RC-aware scoring as the search.
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
                    "score": float(rc_aware_pair_score(df, a, b)),
                }
            )
    rows.sort(key=lambda r: r["score"], reverse=True)
    return rows


def choose_overhangs_by_position(
    df: pd.DataFrame,
    candidates_by_pos: Dict[int, Sequence[DNA]],
    filters: Optional[OverhangFilters] = None,
    beam_width: int = 200,
    *,
    fixed_overhangs: Optional[Sequence[DNA]] = None,
    include_fixed_fixed_in_worst: bool = True,
) -> Tuple[Dict[int, DNA], float]:
    """
    Choose one overhang per position from candidates_by_pos, minimizing the maximum RC-aware cross-talk score
    among all chosen overhangs.

    IMPORTANT (fix for your reported behavior):
      - If fixed_overhangs is provided (e.g., terminal overhangs AATG and GCTT),
        those are included in the minimax objective during optimization.
      - This prevents selecting internal overhangs that massively cross-talk with fixed terminals
        (e.g. AATG vs CATT = 3967) unless literally unavoidable.

    Enforces:
      - chosen overhangs are distinct
      - chosen overhangs cannot be reverse-complements of each other
      - fixed_overhangs must also obey these constraints (with candidates)

    Returns:
      - best_assign: {pos: overhang} for positions in candidates_by_pos
      - best_worst: worst (maximum) pair score among (fixed + chosen),
                   unless include_fixed_fixed_in_worst=False, in which case fixed-fixed
                   pairs do not contribute to best_worst (useful for terminal↔terminal).
    """
    flt = filters or OverhangFilters()

    # Prepare fixed set (dedup, validate)
    fixed_list: List[DNA] = []
    if fixed_overhangs:
        for x in fixed_overhangs:
            x = str(x).upper().strip()
            if len(x) != 4 or any(b not in "ACGT" for b in x):
                raise ValueError(f"Invalid fixed overhang '{x}'. Must be 4 bp A/C/G/T.")
            fixed_list.append(x)

    # Check fixed uniqueness / RC uniqueness
    fixed_set = set(fixed_list)
    if len(fixed_set) != len(fixed_list):
        raise ValueError(f"Duplicate fixed overhangs provided: {fixed_list}")
    fixed_rc_set = set(revcomp(x) for x in fixed_list)
    if fixed_set & fixed_rc_set:
        # This means some fixed oh equals revcomp of another fixed oh (or itself).
        raise ValueError(f"Fixed overhangs contain reverse-complement collision: {fixed_list}")

    # Filter candidates per position
    filtered: Dict[int, List[DNA]] = {}
    for pos, cands in candidates_by_pos.items():
        filtered[pos] = filter_overhangs(cands, flt)
        if not filtered[pos]:
            raise ValueError(f"No feasible overhang candidates for position {pos} after filtering.")

    positions = sorted(filtered.keys())

    def pair(a: DNA, b: DNA) -> float:
        # Hard disallow identical or RC-equal selections (infinite penalty)
        if a == b:
            return float("inf")
        if a == revcomp(b) or revcomp(a) == b:
            return float("inf")
        return float(rc_aware_pair_score(df, a, b))

    # Initial beam includes fixed overhangs already "chosen"
    initial_worst = 0.0
    if include_fixed_fixed_in_worst and len(fixed_list) >= 2:
        for i in range(len(fixed_list)):
            for j in range(i + 1, len(fixed_list)):
                initial_worst = max(initial_worst, pair(fixed_list[i], fixed_list[j]))

    # Track: (assign, chosen, worst_noncrosstalk, min_desired_pair)
    initial_min_desired = float("inf")
    for x in fixed_list:
        initial_min_desired = min(initial_min_desired, desired_pair_score(df, x))

    beam: List[Tuple[Dict[int, DNA], List[DNA], float, float]] = [
        ({}, list(fixed_list), float(initial_worst), float(initial_min_desired))
    ]

    for pos in positions:
        new_beam = []
        for assign, chosen, worst, min_desired in beam:
            chosen_set = set(chosen)
            chosen_rc_set = set(revcomp(x) for x in chosen)

            for cand in candidates_by_pos[pos]:
                cand = str(cand).upper()

                if cand in chosen_set:
                    continue
                if cand in chosen_rc_set:
                    continue

                new_assign = dict(assign)
                new_assign[pos] = cand

                # update worst nonpair crosstalk
                new_worst = worst
                for x in chosen:
                    s = pair(x, cand)
                    if s > new_worst:
                        new_worst = s

                # update "weakest intended pair strength" tie-breaker
                new_min_desired = min(min_desired, desired_pair_score(df, cand))

                new_beam.append((new_assign, chosen + [cand], new_worst, new_min_desired))

        # Objective:
        #  1) minimize worst RC-aware cross-talk among nonpairs (x[2])
        #  2) tie-break: maximize the weakest intended pairing strength (x[3])
        new_beam.sort(key=lambda x: (x[2], -x[3]))
        beam = new_beam[: max(1, int(beam_width))]

        if not beam:
            raise RuntimeError(f"Overhang search failed at position {pos} (beam exhausted).")

    best_assign, chosen_all, best_worst, best_min_desired = beam[0]
    return best_assign, float(best_worst)