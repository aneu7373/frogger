
from __future__ import annotations

from typing import Any, Mapping


def _w(weights: Mapping[str, float], key: str, default: float = 0.0) -> float:
    return float(weights.get(key, default)) if weights is not None else default


def _gc_bounds_penalty(gc: float, gc_min: float | None, gc_max: float | None) -> float:
    """
    Penalty is the magnitude of violation outside [gc_min, gc_max]; 0 if within bounds.
    """
    penalty = 0.0
    if gc_min is not None and gc < gc_min:
        penalty += (gc_min - gc)
    if gc_max is not None and gc > gc_max:
        penalty += (gc - gc_max)
    return penalty


def score_candidate(metrics: Mapping[str, Any], weights: Mapping[str, float]) -> float:
    """
    Composite score (higher is better).

    Metrics expected:
      gc, mean_codon_freq, entropy, rna_5p, rna_full, restriction_hits, motif_hits,
      forbidden_hits (optional)

    Weights:
      gc, gc_target, gc_min, gc_max, codon, entropy, rna_5p, rna_full, restriction, motif, forbidden
    """
    score = 0.0

    gc = float(metrics.get("gc", 0.0))
    gc_min = weights.get("gc_min")
    gc_max = weights.get("gc_max")
    if gc_min is not None or gc_max is not None:
        gc_pen = _gc_bounds_penalty(gc, gc_min, gc_max)
        score -= _w(weights, "gc", 1.0) * gc_pen
    else:
        gc_tgt = float(weights.get("gc_target", 0.5))
        score -= _w(weights, "gc", 1.0) * abs(gc - gc_tgt)

    score += _w(weights, "codon", 1.0) * float(metrics.get("mean_codon_freq", 0.0))
    score += _w(weights, "entropy", 0.0) * float(metrics.get("entropy", 0.0))

    # RNA terms: sign convention controlled by weights
    score += _w(weights, "rna_5p", 0.0) * float(metrics.get("rna_5p", 0.0))
    score += _w(weights, "rna_full", 0.0) * float(metrics.get("rna_full", 0.0))

    score -= _w(weights, "restriction", 0.0) * float(metrics.get("restriction_hits", 0))
    score -= _w(weights, "motif", 0.0) * float(metrics.get("motif_hits", 0))
    score -= _w(weights, "forbidden", 0.0) * float(metrics.get("forbidden_hits", 0))

    return score

