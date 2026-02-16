# codonopt/core/optimizer.py

from __future__ import annotations

import random
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional

from codonopt.core.generation import filter_codons_for_constraints
from codonopt.exceptions import ConstraintError
from codonopt.utils import calculate_gc


def _gc_count(seq: str) -> int:
    """Count G/C bases in a DNA string (case-insensitive)."""
    if not seq:
        return 0
    s = seq.upper()
    return s.count("G") + s.count("C")


def _is_gc_feasible(
    gc_so_far: int,
    next_pos: int,
    total_nt: int,
    suffix_min_gc: List[int],
    suffix_max_gc: List[int],
    gc_min: Optional[float],
    gc_max: Optional[float],
    eps: float = 1e-6,
) -> bool:
    """
    Feasibility-based GC pruning.

    Only reject a prefix if it is mathematically impossible to end within [gc_min, gc_max]
    given per-position min/max GC achievable from remaining amino acids (optimistic bounds).
    """
    if gc_min is None and gc_max is None:
        return True

    min_total_gc = gc_so_far + suffix_min_gc[next_pos]
    max_total_gc = gc_so_far + suffix_max_gc[next_pos]

    min_frac = min_total_gc / float(total_nt)
    max_frac = max_total_gc / float(total_nt)

    if gc_min is not None and max_frac < (gc_min - eps):
        return False
    if gc_max is not None and min_frac > (gc_max + eps):
        return False
    return True


def _weighted_choice(rng: random.Random, codons: List[str], weights: List[float]) -> str:
    if not codons:
        raise ConstraintError("No candidate codons available after filtering.")
    if len(codons) != len(weights):
        raise ConstraintError("Internal error: codons/weights length mismatch.")
    return rng.choices(codons, weights=weights, k=1)[0]


def _weighted_permutation(rng: random.Random, codons: List[str], weights: List[float]) -> List[str]:
    """
    Weighted random ordering (without replacement), biased to try higher-weight codons earlier.

    Efraimidis-Spirakis style keys: key = U^(1/w)
    Higher weights tend to yield larger keys, thus appear earlier when sorting descending.
    """
    if not codons:
        return []
    if len(codons) != len(weights):
        raise ConstraintError("Internal error: codons/weights length mismatch.")

    items = []
    for c, w in zip(codons, weights):
        w = float(w)
        if w <= 0:
            key = -1.0
        else:
            u = rng.random()
            key = u ** (1.0 / w)
        items.append((key, c))

    items.sort(reverse=True, key=lambda x: x[0])
    return [c for _, c in items]


@dataclass
class OptimizerDebug:
    """Lightweight debug counters; included in failure messages."""
    # filtering / constraint related
    filter_calls: int = 0
    filter_failures: int = 0  # filter raised "no valid codons"
    # gc related
    gc_pruned_infeasible: int = 0
    gc_pruned_final_out_of_range: int = 0
    # search mechanics
    strict_position_retries: int = 0
    kleinbub_backtracks: int = 0
    kleinbub_steps: int = 0
    # misc
    notes: List[str] = field(default_factory=list)

    def summary(self) -> str:
        return (
            "OptimizerDebug("
            f"filter_calls={self.filter_calls}, "
            f"filter_failures={self.filter_failures}, "
            f"gc_pruned_infeasible={self.gc_pruned_infeasible}, "
            f"gc_pruned_final_out_of_range={self.gc_pruned_final_out_of_range}, "
            f"strict_position_retries={self.strict_position_retries}, "
            f"kleinbub_backtracks={self.kleinbub_backtracks}, "
            f"kleinbub_steps={self.kleinbub_steps}"
            ")"
        )


def optimize_sequence(
    dna=None,
    protein=None,
    codon_table=None,
    avoid_codons=None,
    avoid_motifs=None,
    max_homopolymer=5,
    gc_min=None,
    gc_max=None,
    logger=None,
    optimization_mode="kleinbub",
    max_attempts=1000,
    backtrack_window=10,
    min_codon_fraction=0.05,
    seed: Optional[int] = None,
    debug: Optional[Dict] = None,
):
    """
    Deterministic codon optimization (when seed is provided).

    optimization_mode:
      - 'strict': greedy, retries per position, fast
      - 'kleinbub': bounded backtracking search, higher yield, slower

    codon_table must be:
      dict: AA -> list of (codon, fraction)

    min_codon_fraction:
      excludes codons for an AA if fraction < cutoff.
      default 0.05.

    debug:
      if provided (dict-like), will be populated with debug counters under key 'optimizer'.
    """
    rng = random.Random(seed) if seed is not None else random.Random()

    dbg = OptimizerDebug()
    if debug is not None:
        debug["optimizer"] = dbg  # caller can read dbg.summary() or fields

    avoid_codons = avoid_codons or []
    avoid_motifs = avoid_motifs or []

    if protein is None:
        raise ConstraintError("optimize_sequence requires protein input (DNA is translated upstream).")
    if codon_table is None:
        raise ConstraintError("codon_table is required.")

    mode = (optimization_mode or "kleinbub").strip().lower()
    if mode not in ("strict", "kleinbub"):
        raise ConstraintError("optimization_mode must be 'strict' or 'kleinbub'")

    try:
        min_cf = float(min_codon_fraction)
    except Exception:
        raise ConstraintError(f"min_codon_fraction must be a float; got {min_codon_fraction!r}")
    if min_cf < 0 or min_cf > 1:
        raise ConstraintError(f"min_codon_fraction must be between 0 and 1; got {min_cf}")

    # Normalize/validate GC bounds
    if gc_min is not None:
        try:
            gc_min = float(gc_min)
        except Exception:
            raise ConstraintError(f"gc_min must be a float in [0,1]; got {gc_min!r}")
        if not (0.0 <= gc_min <= 1.0):
            raise ConstraintError(f"gc_min must be between 0 and 1; got {gc_min}")

    if gc_max is not None:
        try:
            gc_max = float(gc_max)
        except Exception:
            raise ConstraintError(f"gc_max must be a float in [0,1]; got {gc_max!r}")
        if not (0.0 <= gc_max <= 1.0):
            raise ConstraintError(f"gc_max must be between 0 and 1; got {gc_max}")

    if (gc_min is not None) and (gc_max is not None) and (gc_min > gc_max):
        raise ConstraintError(f"gc_min ({gc_min}) cannot be greater than gc_max ({gc_max}).")

    # Build AA -> (codons, weights) after rarity cutoff
    aa_codons: Dict[str, Tuple[List[str], List[float]]] = {}
    for aa in set(protein):
        rows = codon_table.get(aa, [])
        if not rows:
            raise ConstraintError(f"No codons available in codon table for amino acid '{aa}'")
        filtered = [(c, float(fr)) for c, fr in rows if float(fr) >= min_cf]
        if not filtered:
            raise ConstraintError(
                f"All codons for amino acid '{aa}' are below min_codon_fraction={min_cf}. "
                f"Lower the cutoff or change the codon table."
            )
        codons = [c for c, _ in filtered]
        weights = [fr for _, fr in filtered]
        aa_codons[aa] = (codons, weights)

    # Precompute optimistic per-position GC bounds for feasibility pruning
    total_nt = 3 * len(protein)
    min_gc_per_aa: Dict[str, int] = {}
    max_gc_per_aa: Dict[str, int] = {}
    for aa, (codons, _weights) in aa_codons.items():
        gc_counts = [_gc_count(c) for c in codons]
        min_gc_per_aa[aa] = min(gc_counts)
        max_gc_per_aa[aa] = max(gc_counts)

    suffix_min_gc = [0] * (len(protein) + 1)
    suffix_max_gc = [0] * (len(protein) + 1)
    for idx in range(len(protein) - 1, -1, -1):
        aa = protein[idx]
        suffix_min_gc[idx] = suffix_min_gc[idx + 1] + min_gc_per_aa[aa]
        suffix_max_gc[idx] = suffix_max_gc[idx + 1] + max_gc_per_aa[aa]

    eps = 1e-6

    # ---------------------------
    # STRICT MODE
    # ---------------------------
    if mode == "strict":
        seq = ""
        for pos, aa in enumerate(protein):
            codons, weights = aa_codons[aa]
            attempts = 0

            while attempts < max_attempts:
                dbg.filter_calls += 1
                try:
                    valid_codons = filter_codons_for_constraints(
                        codons=codons,
                        current_seq=seq,
                        avoid_codons=avoid_codons,
                        avoid_motifs=avoid_motifs,
                        max_homopolymer=max_homopolymer,
                        context=f"strict pos={pos} aa={aa}",
                    )
                except ConstraintError as e:
                    dbg.filter_failures += 1
                    raise ConstraintError(
                        f"Strict optimization failed: no valid codons at pos={pos} aa={aa}. "
                        f"Prefix_len_nt={len(seq)}. Underlying: {e}. {dbg.summary()}"
                    ) from e

                wmap = {c: w for c, w in zip(codons, weights)}
                valid_weights = [wmap[c] for c in valid_codons]
                codon = _weighted_choice(rng, valid_codons, valid_weights)
                candidate = seq + codon

                cand_gc = _gc_count(candidate)
                is_last = (pos == len(protein) - 1)

                if is_last:
                    gc = calculate_gc(candidate)
                    if (gc_min is not None and gc < (gc_min - eps)) or (gc_max is not None and gc > (gc_max + eps)):
                        dbg.gc_pruned_final_out_of_range += 1
                        attempts += 1
                        dbg.strict_position_retries += 1
                        continue
                    seq = candidate
                    if logger:
                        logger.debug(f"{aa} -> {codon} (strict)")
                    break
                else:
                    if not _is_gc_feasible(
                        gc_so_far=cand_gc,
                        next_pos=pos + 1,
                        total_nt=total_nt,
                        suffix_min_gc=suffix_min_gc,
                        suffix_max_gc=suffix_max_gc,
                        gc_min=gc_min,
                        gc_max=gc_max,
                        eps=eps,
                    ):
                        dbg.gc_pruned_infeasible += 1
                        attempts += 1
                        dbg.strict_position_retries += 1
                        continue

                    seq = candidate
                    if logger:
                        logger.debug(f"{aa} -> {codon} (strict)")
                    break

            else:
                raise ConstraintError(
                    f"Strict optimization failed after {max_attempts} attempts at pos={pos} aa={aa}. "
                    f"Prefix_len_nt={len(seq)}. {dbg.summary()}"
                )

        # Final guard (redundant but safe)
        final_gc = calculate_gc(seq)
        if (gc_min is not None and final_gc < (gc_min - eps)) or (gc_max is not None and final_gc > (gc_max + eps)):
            raise ConstraintError(
                f"Strict optimization produced sequence outside GC bounds unexpectedly "
                f"(gc={final_gc}, bounds=[{gc_min},{gc_max}]). {dbg.summary()}"
            )
        return seq

    # ---------------------------
    # KLEINBUB MODE (bounded backtracking)
    # ---------------------------
    if backtrack_window < 1:
        backtrack_window = 1

    chosen: List[Tuple[str, str]] = []  # (aa, codon)
    seq = ""
    i = 0
    total_steps = 0

    search_limit = int(max_attempts)

    candidate_lists: List[Optional[List[str]]] = [None] * len(protein)
    candidate_index = [0] * len(protein)

    while i < len(protein):
        total_steps += 1
        dbg.kleinbub_steps = total_steps

        if total_steps > search_limit:
            raise ConstraintError(
                f"Kleinbub optimization exceeded search limit ({search_limit}). "
                f"Prefix_len_nt={len(seq)}. {dbg.summary()}"
            )

        aa = protein[i]
        codons, weights = aa_codons[aa]

        if candidate_lists[i] is None:
            dbg.filter_calls += 1
            try:
                valid_codons = filter_codons_for_constraints(
                    codons=codons,
                    current_seq=seq,
                    avoid_codons=avoid_codons,
                    avoid_motifs=avoid_motifs,
                    max_homopolymer=max_homopolymer,
                    context=f"kleinbub pos={i} aa={aa}",
                )
            except ConstraintError:
                dbg.filter_failures += 1
                candidate_lists[i] = []
                candidate_index[i] = 0
            else:
                wmap = {c: w for c, w in zip(codons, weights)}
                valid_weights = [wmap[c] for c in valid_codons]
                ordered = _weighted_permutation(rng, valid_codons, valid_weights)
                candidate_lists[i] = ordered
                candidate_index[i] = 0

        # exhausted candidates -> backtrack
        if candidate_index[i] >= len(candidate_lists[i]):
            candidate_lists[i] = None
            candidate_index[i] = 0

            if i == 0:
                raise ConstraintError(
                    "Kleinbub optimization could not satisfy constraints from the start. "
                    f"{dbg.summary()}"
                )

            back_to = max(0, i - backtrack_window)
            dbg.kleinbub_backtracks += 1

            while len(chosen) > back_to:
                chosen.pop()

            seq = "".join(c for (_, c) in chosen)

            for j in range(back_to, len(protein)):
                candidate_lists[j] = None
                candidate_index[j] = 0

            i = back_to
            continue

        codon = candidate_lists[i][candidate_index[i]]
        candidate_index[i] += 1

        candidate_seq = seq + codon
        cand_gc = _gc_count(candidate_seq)
        is_last = (i == len(protein) - 1)

        if is_last:
            gc = calculate_gc(candidate_seq)
            if (gc_min is not None and gc < (gc_min - eps)) or (gc_max is not None and gc > (gc_max + eps)):
                dbg.gc_pruned_final_out_of_range += 1
                continue
        else:
            if not _is_gc_feasible(
                gc_so_far=cand_gc,
                next_pos=i + 1,
                total_nt=total_nt,
                suffix_min_gc=suffix_min_gc,
                suffix_max_gc=suffix_max_gc,
                gc_min=gc_min,
                gc_max=gc_max,
                eps=eps,
            ):
                dbg.gc_pruned_infeasible += 1
                continue

        chosen.append((aa, codon))
        seq = candidate_seq
        if logger:
            logger.debug(f"{aa} -> {codon} (kleinbub)")
        i += 1

    final_gc = calculate_gc(seq)
    if (gc_min is not None and final_gc < (gc_min - eps)) or (gc_max is not None and final_gc > (gc_max + eps)):
        raise ConstraintError(
            f"Kleinbub optimization produced sequence outside GC bounds unexpectedly "
            f"(gc={final_gc}, bounds=[{gc_min},{gc_max}]). {dbg.summary()}"
        )
    return seq
