from __future__ import annotations

import random
from pathlib import Path
from typing import Dict, List, Set, Optional

DNA = str


def revcomp(seq: DNA) -> DNA:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def gc_fraction(seq: DNA) -> float:
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0


def max_homopolymer_run(seq: DNA) -> int:
    if not seq:
        return 0
    best = cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
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
            out.add(s[i : i + k])
        rc = revcomp(s)
        for i in range(0, len(rc) - k + 1):
            out.add(rc[i : i + k])
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

    Ensures primers do not contain any k-mer present in avoid_kmers
    (forward/rc handled already in build_avoid_kmers).
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
        for i in range(0, len(pr) - k + 1):
            if pr[i : i + k] in avoid_kmers:
                return False
        if pr in used or revcomp(pr) in used:
            return False
        return True

    for pos in range(1, n_positions + 1):
        pair: Dict[str, str] = {}
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


def load_barcodes_fasta(path: Path) -> Dict[str, DNA]:
    """
    Load barcode sequences from a FASTA file.
    Uses the first token of the FASTA header (after '>') as the barcode ID.

    Example:
      >901
      ACTG...
      >902 some description
      TGCA...

    Returns { "901": "ACTG...", "902": "TGCA...", ... } (uppercased, A/C/G/T only).
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Barcode FASTA not found: {p}")
    out: Dict[str, DNA] = {}
    cur_id: Optional[str] = None
    cur_seq: List[str] = []
    with p.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seq = "".join(cur_seq).upper().replace("U", "T")
                    seq = "".join([c for c in seq if c in "ACGT"])
                    if seq:
                        out[cur_id] = seq
                cur_id = line[1:].strip().split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id is not None:
            seq = "".join(cur_seq).upper().replace("U", "T")
            seq = "".join([c for c in seq if c in "ACGT"])
            if seq:
                out[cur_id] = seq
    return out


def primers_from_barcodes(
    n_positions: int,
    barcode_fasta: Path,
    forward_id: str = "901",
    reverse_id: str = "902",
    forward_seq: Optional[str] = None,
    reverse_seq: Optional[str] = None,
    reverse_is_revcomp: bool = True,
) -> Dict[int, Dict[str, str]]:
    """
    Build per-position primers from barcode IDs (or explicit sequences).

    Default behavior:
      - left primer  = barcode 901 (forward)
      - right primer = reverse-complement(barcode 902)

    If forward_seq/reverse_seq are provided, those override the FASTA lookup.
    """
    bcs = load_barcodes_fasta(Path(barcode_fasta))

    if forward_seq is None:
        if str(forward_id) not in bcs:
            raise ValueError(f"Forward barcode id '{forward_id}' not found in {barcode_fasta}")
        forward_seq = bcs[str(forward_id)]
    if reverse_seq is None:
        if str(reverse_id) not in bcs:
            raise ValueError(f"Reverse barcode id '{reverse_id}' not found in {barcode_fasta}")
        reverse_seq = bcs[str(reverse_id)]

    forward_seq = str(forward_seq).upper().replace("U", "T")
    reverse_seq = str(reverse_seq).upper().replace("U", "T")
    if reverse_is_revcomp:
        reverse_seq = revcomp(reverse_seq)

    out: Dict[int, Dict[str, str]] = {}
    for pos in range(1, n_positions + 1):
        out[pos] = {"left": forward_seq, "right": reverse_seq}
    return out
