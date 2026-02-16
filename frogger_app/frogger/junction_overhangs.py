from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from frogger.common.codon_table import CodonEntry


DNA = str


@dataclass(frozen=True)
class JunctionChoice:
    """
    Specifies a concrete codon pair and a 4-mer window within the 6 nt across the boundary.
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
) -> Dict[str, JunctionChoice]:
    """
    For a boundary between aa_left and aa_right, enumerate all possible 4-mer overhangs
    that can appear within the 6 nt formed by (codon_left + codon_right), considering synonymous codons.

    Returns mapping overhang -> best JunctionChoice (best by score).
    """
    avoid = set([c.upper() for c in (avoid_codons or [])])
    if aa_left not in codons or aa_right not in codons:
        raise ValueError(f"Missing codons for AA boundary {aa_left}-{aa_right} in codon table.")

    best: Dict[str, JunctionChoice] = {}
    for cl in codons[aa_left]:
        if cl.codon in avoid:
            continue
        for cr in codons[aa_right]:
            if cr.codon in avoid:
                continue
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
