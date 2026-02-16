import re
from codonopt.exceptions import ConstraintError


def violates_homopolymer(seq, max_run):
    runs = re.findall(r"(A+|T+|G+|C+)", seq)
    return any(len(r) > max_run for r in runs)


def filter_codons_for_constraints(
    codons,
    current_seq,
    avoid_codons,
    avoid_motifs,
    max_homopolymer,
    context: str | None = None,
):
    """
    Filter candidate codons using local constraints:
      - max homopolymer (lookback window)
      - avoid_codons
      - avoid_motifs (substring search across current_seq+codon)

    Raises a ConstraintError with informative rejection counts if all codons are filtered out.
    """
    rejected_homopolymer = 0
    rejected_avoid_codon = 0
    rejected_avoid_motif = 0

    valid = []
    # only need the tail for homopolymer checking
    tail = current_seq[-(max_homopolymer - 1) :] if max_homopolymer and max_homopolymer > 1 else ""

    for codon in codons:
        test_seq = tail + codon
        if violates_homopolymer(test_seq, max_homopolymer):
            rejected_homopolymer += 1
            continue
        if codon in avoid_codons:
            rejected_avoid_codon += 1
            continue
        # motif check over the full concatenation
        joined = current_seq + codon
        if any(m in joined for m in avoid_motifs):
            rejected_avoid_motif += 1
            continue
        valid.append(codon)

    if not valid:
        ctx = f" ({context})" if context else ""
        raise ConstraintError(
            "No valid codons remaining due to constraints"
            f"{ctx}. Rejections: homopolymer={rejected_homopolymer}, "
            f"avoid_codons={rejected_avoid_codon}, avoid_motifs={rejected_avoid_motif}. "
            f"Prefix_tail='{tail}'. Total_candidates={len(codons)}"
        )

    return valid
