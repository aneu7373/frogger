
from __future__ import annotations

from typing import Optional
from Bio.Seq import Seq


def translate_dna(dna: str, *, to_stop: bool = False, table: Optional[int] = None) -> str:
    """
    Translate a DNA coding sequence to amino acids (one-letter). May include '*'
    if to_stop=False and internal stops are present.
    """
    seq = Seq(dna.upper())
    if table is None:
        aa = seq.translate(to_stop=to_stop)
    else:
        aa = seq.translate(table=table, to_stop=to_stop)
    return str(aa)


def validate_orf(candidate_dna: str, reference_aa: str, *, table: Optional[int] = None) -> bool:
    """
    True if translating candidate_dna yields exactly reference_aa.
    """
    aa = translate_dna(candidate_dna, to_stop=False, table=table)
    return aa == reference_aa

