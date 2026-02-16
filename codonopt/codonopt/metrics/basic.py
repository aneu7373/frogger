
from __future__ import annotations

import math
from typing import Iterable


def gc_content(seq: str) -> float:
    """
    GC fraction in [0,1]. Case-insensitive; returns 0.0 for empty input.
    """
    if not seq:
        return 0.0
    s = seq.upper()
    gc = s.count("G") + s.count("C")
    return gc / len(s)


def shannon_entropy(values: Iterable[float]) -> float:
    """
    Shannon entropy H = -sum_i p_i log2 p_i.
    If inputs are not normalized, they are normalized internally.
    """
    vals = list(values)
    total = sum(vals)
    if total <= 0:
        return 0.0
    entropy = 0.0
    for v in vals:
        p = v / total
        if p > 0.0:
            entropy -= p * math.log2(p)
    return entropy

