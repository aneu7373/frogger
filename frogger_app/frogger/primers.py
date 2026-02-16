from __future__ import annotations
import random
from typing import Dict, List, Set, Tuple

DNA = str

def revcomp(seq: DNA) -> DNA:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]

def gc_fraction(seq: DNA) -> float:
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s)

def max_homopolymer_run(seq: DNA) -> int:
    best = cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best

def build_avoid_kmers(seqs: List[DNA], k: int) -> Set[DNA]:
    out: Set[DNA] = set()
    for s in seqs:
        s = s.upper()
        for i in range(0, len(s) - k + 1):
            out.add(s[i:i+k])
        rc = revcomp(s)
        for i in range(0, len(rc) - k + 1):
            out.add(rc[i:i+k])
    return out

def _rand_dna(rng: random.Random, n: int) -> DNA:
    return "".join(rng.choice("ACGT") for _ in range(n))

def generate_primers_by_pos(
    n_positions: int,
    length: int,
    avoid_kmers: Set[DNA],
    k: int,
    min_gc: float = 0.35,
    max_gc: float = 0.65,
    max_homopolymer: int = 3,
    avoid_motifs: List[DNA] | None = None,
    seed: int = 1337,
    max_tries: int = 200000,
) -> Dict[int, Dict[str, str]]:
    """
    Returns {pos: {left: primer, right: primer}}.
    Ensures primers do not contain any k-mer present in avoid_kmers (forward/rc handled already).
    """
    rng = random.Random(seed)
    avoid_motifs = [m.upper() for m in (avoid_motifs or [])]

    used: Set[DNA] = set()
    out: Dict[int, Dict[str, str]] = {}

    def ok(pr: DNA) -> bool:
        pr = pr.upper()
        if not (min_gc <= gc_fraction(pr) <= max_gc):
            return False
        if max_homopolymer_run(pr) > max_homopolymer:
            return False
        for m in avoid_motifs:
            if m and m in pr:
                return False
        # screen k-mers inside primer
        for i in range(0, len(pr) - k + 1):
            if pr[i:i+k] in avoid_kmers:
                return False
        # also prevent exact reuse / reverse complement reuse
        if pr in used or revcomp(pr) in used:
            return False
        return True

    for pos in range(1, n_positions + 1):
        pair = {}
        for side in ("left", "right"):
            found = None
            tries = 0
            while tries < max_tries:
                tries += 1
                cand = _rand_dna(rng, length)
                if ok(cand):
                    found = cand
                    used.add(found)
                    break
            if found is None:
                raise RuntimeError(f"Failed to generate primer for pos={pos} side={side} after {max_tries} tries.")
            pair[side] = found
        out[pos] = pair

    return out
