from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import pandas as pd


@dataclass(frozen=True)
class CodonEntry:
    codon: str
    fraction: float


def load_codon_usage_xlsx(path: Path, sheet_name: Optional[str] = None) -> Dict[str, List[CodonEntry]]:
    """
    Load a simple codon-usage table with columns: AA, Codon, Fraction.
    Returns mapping AA(one-letter) -> list of CodonEntry, sorted by descending Fraction.
    Stops ('*') are excluded.
    """
    xls = pd.ExcelFile(path)
    sheet = sheet_name or xls.sheet_names[0]
    df = pd.read_excel(path, sheet_name=sheet)
    required = {"AA", "Codon", "Fraction"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Codon table missing columns: {sorted(missing)}")
    out: Dict[str, List[CodonEntry]] = {}
    for _, row in df.iterrows():
        aa = str(row["AA"]).strip()
        codon = str(row["Codon"]).strip().upper()
        frac = float(row["Fraction"])
        if aa == "*" or aa == "STOP":
            continue
        if len(codon) != 3 or any(b not in "ACGT" for b in codon):
            continue
        out.setdefault(aa, []).append(CodonEntry(codon=codon, fraction=frac))
    for aa, lst in out.items():
        out[aa] = sorted(lst, key=lambda x: x.fraction, reverse=True)
    return out
