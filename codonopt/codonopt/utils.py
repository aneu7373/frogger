from __future__ import annotations

import re
from typing import Set, Tuple, Optional

from Bio.Seq import Seq

from codonopt.exceptions import InputFormatError

DNA_CHARS: Set[str] = set("ACGTUN")  # allow U (RNA) and N
AA_CHARS: Set[str] = set("ACDEFGHIKLMNPQRSTVWY")  # 20 canonical AAs


def is_dna_sequence(seq: str) -> bool:
    """
    Strict-ish DNA/RNA detector:
      - if every character is in A/C/G/T/U/N (case-insensitive), treat as nucleotide
      - otherwise treat as protein-ish (or invalid input)
    """
    s = (seq or "").strip().upper()
    if not s:
        return False
    return all(c in DNA_CHARS for c in s)


def is_protein_sequence(seq: str) -> bool:
    """
    Protein detector:
      - returns True if sequence is NOT purely nucleotide alphabet (ACGTUN)
    """
    return not is_dna_sequence(seq)


def calculate_gc(seq: str) -> float:
    seq = (seq or "").upper()
    if not seq:
        return 0.0
    gc_count = seq.count("G") + seq.count("C")
    return gc_count / len(seq)


def max_homopolymer_length(seq: str) -> int:
    seq = (seq or "").upper()
    runs = re.findall(r"(A+|T+|G+|C+)", seq)
    return max((len(r) for r in runs), default=0)


def translate_cds_to_protein(dna: str) -> str:
    """
    Translate a CDS using the standard genetic code.

    - Accepts DNA or RNA characters; converts U->T
    - Requires length divisible by 3
    - Translates full length (to_stop=False)
    - If the result ends with '*', strips the terminal stop for comparison purposes
      (because codonopt optimizes the coding region only).
    """
    if dna is None:
        raise InputFormatError("translate_cds_to_protein: dna is None")

    s = str(dna).strip().upper().replace("U", "T")
    if not s:
        raise InputFormatError("translate_cds_to_protein: empty DNA sequence")

    if len(s) % 3 != 0:
        raise InputFormatError(f"CDS length not divisible by 3 (len={len(s)})")

    aa = str(Seq(s).translate(to_stop=False))

    # Strip terminal stop for comparison
    if aa.endswith("*"):
        aa = aa[:-1]

    return aa


def verify_backtranslation_matches_protein(
    candidate_dna: str,
    target_protein: str,
    *,
    allow_terminal_stop: bool = True,
) -> Tuple[bool, str, Optional[str]]:
    """
    Verify that candidate_dna translates back to exactly target_protein.

    Returns:
      (ok, reason, translated_aa)

    Behavior:
      - Rejects if length % 3 != 0
      - Rejects if translated AA contains internal stop codons
      - Strips terminal '*' if allow_terminal_stop=True
      - Requires exact string equality with target_protein
    """
    try:
        aa = translate_cds_to_protein(candidate_dna)
    except Exception as e:
        return False, f"translation_failed: {e}", None

    # If BioPython translation produced internal stops, they'll still be present
    # (terminal stop stripped in translate_cds_to_protein).
    if "*" in aa:
        return False, "internal_stop_codon_in_candidate", aa

    tgt = (target_protein or "").strip().upper()
    if not tgt:
        return False, "empty_target_protein", aa

    if aa != tgt:
        # Find first mismatch index (helps debugging)
        mismatch_at = None
        for i, (a, b) in enumerate(zip(aa, tgt)):
            if a != b:
                mismatch_at = i
                break
        if mismatch_at is None and len(aa) != len(tgt):
            mismatch_at = min(len(aa), len(tgt))

        return False, f"aa_mismatch_at:{mismatch_at}", aa

    return True, "", aa
