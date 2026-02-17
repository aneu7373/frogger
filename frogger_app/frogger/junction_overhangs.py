from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from frogger.common.codon_table import CodonEntry


DNA = str


@dataclass(frozen=True)
class JunctionChoice:
    """
    Specifies a concrete codon pair and a 4-mer window within the 6 nt across the boundary.

    NOTE: In unified FROGGER, junction overhang enumeration is driven by the amino-acid
    context *downstream* of the cut (see packager.select_global_overhangs_and_forced_codons).
    This helper simply enumerates feasible 4-mers from two adjacent amino acids given a
    codon-usage table.
    """
    left_codon: str
    right_codon: str
    overhang: str
    window_offset: int  # 0,1,2 representing 6nt[offset:offset+4]
    score: float        # higher is better (e.g., codon usage product)


def _windows_4mer(six_nt: str) -> List[Tuple[int, str]]:
    return [(0, six_nt[0:4]), (1, six_nt[1:5]), (2, six_nt[2:6])]


def enumerate_boundary_overhangs(
    aa_left: str,
    aa_right: str,
    codons: Dict[str, List[CodonEntry]],
    avoid_codons: Optional[Sequence[str]] = None,
    top_k_codons: int = 2,
) -> Dict[str, JunctionChoice]:
    """
    For a boundary between aa_left and aa_right, enumerate all possible 4-mer overhangs
    that can appear within the 6 nt formed by (codon_left + codon_right), considering
    synonymous codons.

    Key update:
      - Restrict synonymous-codon enumeration to the top-K most frequent codons per AA
        (default K=2). This keeps the overhang search focused on highly used codons.

    Returns mapping: overhang -> best JunctionChoice (best by score).
    """
    avoid = set([c.upper() for c in (avoid_codons or [])])
    if aa_left not in codons or aa_right not in codons:
        raise ValueError(f"Missing codons for AA boundary {aa_left}-{aa_right} in codon table.")

    k = max(1, int(top_k_codons))
    left_list = [ce for ce in codons[aa_left] if ce.codon.upper() not in avoid][:k]
    right_list = [ce for ce in codons[aa_right] if ce.codon.upper() not in avoid][:k]

    if not left_list or not right_list:
        raise ValueError(
            f"No usable codons after filtering for AA boundary {aa_left}-{aa_right} "
            f"(avoid_codons may be too strict)."
        )

    best: Dict[str, JunctionChoice] = {}
    for cl in left_list:
        for cr in right_list:
            six = (cl.codon + cr.codon).upper()
            base_score = float(cl.fraction) * float(cr.fraction)
            for off, oh in _windows_4mer(six):
                if any(b not in "ACGT" for b in oh):
                    continue
                prev = best.get(oh)
                if prev is None or base_score > prev.score:
                    best[oh] = JunctionChoice(
                        left_codon=cl.codon,
                        right_codon=cr.codon,
                        overhang=oh,
                        window_offset=off,
                        score=base_score,
                    )
    return best


def intersect_candidate_sets(
    per_gene_candidates: List[Dict[str, JunctionChoice]]
) -> Dict[str, JunctionChoice]:
    """
    Given per-gene mapping overhang -> best choice for that gene, return intersection of overhang keys.
    The returned JunctionChoice is a placeholder (from the first gene); per-gene choices are retained separately.
    """
    if not per_gene_candidates:
        return {}

    keys = set(per_gene_candidates[0].keys())
    for d in per_gene_candidates[1:]:
        keys &= set(d.keys())
    return {k: per_gene_candidates[0][k] for k in sorted(keys)}
