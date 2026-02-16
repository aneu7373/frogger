
from __future__ import annotations

from typing import Iterable


def count_hits(
    seq: str,
    motifs: Iterable[str],
    *,
    overlapping: bool = True,
    case_insensitive: bool = True,
) -> int:
    """
    Count motif occurrences in seq. Overlapping hits allowed by default.

    Parameters
    ----------
    overlapping : bool
        If True, counts overlapping hits (e.g., 'AAA' in 'AAAA' -> 2).
    case_insensitive : bool
        If True, uppercases sequence and motifs before searching.
    """
    if not seq:
        return 0
    s = seq.upper() if case_insensitive else seq
    hits = 0
    for motif in motifs:
        if not motif:
            continue
        m = motif.upper() if case_insensitive else motif
        start = 0
        step = 1 if overlapping else max(1, len(m))
        while True:
            idx = s.find(m, start)
            if idx == -1:
                break
            hits += 1
            start = idx + step
    return hits

