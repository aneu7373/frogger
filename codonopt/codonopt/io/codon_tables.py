# codonopt/io/codon_tables.py

from __future__ import annotations

from typing import Dict, List, Tuple, Optional
import re

import pandas as pd

from codonopt.exceptions import InputFormatError


CodonTable = Dict[str, List[Tuple[str, float]]]


def _clean_colname(c: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(c).strip().lower())


def _norm_aa(x: str) -> str:
    s = str(x).strip().upper()
    if s in ("STOP", "STOPCODON", "TER", "TERM", "*"):
        return "*"
    if len(s) != 1:
        # allow e.g. "Ala" -> take first char only is unsafe; reject
        raise InputFormatError(f"Invalid amino acid value: {x!r} (expected single-letter AA or '*')")
    return s


def _norm_codon(x: str) -> str:
    s = str(x).strip().upper().replace("U", "T")
    if not re.fullmatch(r"[ACGT]{3}", s):
        raise InputFormatError(f"Invalid codon: {x!r} (expected A/C/G/T triplet)")
    return s


def _parse_codons_cell(cell: str) -> List[str]:
    # allow comma, pipe, whitespace separators
    s = str(cell).strip()
    if not s:
        return []
    parts = re.split(r"[,\|\s]+", s)
    codons = []
    for p in parts:
        p = p.strip()
        if not p:
            continue
        codons.append(_norm_codon(p))
    return codons


def _normalize_fractions(rows: List[Tuple[str, float]]) -> List[Tuple[str, float]]:
    total = sum(max(0.0, float(fr)) for _, fr in rows)
    if total <= 0:
        # fallback: uniform if all 0 or missing
        n = len(rows)
        if n == 0:
            return []
        return [(c, 1.0 / n) for c, _ in rows]
    return [(c, float(fr) / total) for c, fr in rows]


def load_codon_table_from_xlsx(path: str, sheet_name: Optional[str] = None) -> CodonTable:
    """
    Load codon table from an XLSX file.

    Returns:
      dict: AA -> list of (codon, fraction), fractions normalized per AA.

    Supports:
      - Long format with columns like:
          AA / Amino Acid, Codon, Fraction / Frequency / Usage
      - Wide format with columns like:
          Amino Acid, Codons
        (fractions assumed uniform within AA)
    """
    try:
        xl = pd.ExcelFile(path, engine="openpyxl")
    except Exception as e:
        raise InputFormatError(f"Failed to open codon table XLSX: {path}. Error: {e}")

    sheets = xl.sheet_names
    if not sheets:
        raise InputFormatError(f"No sheets found in codon table XLSX: {path}")

    if sheet_name is None:
        # pick the first sheet that has any rows
        chosen = None
        for s in sheets:
            df = xl.parse(s)
            if df is not None and not df.empty:
                chosen = s
                break
        if chosen is None:
            raise InputFormatError(f"All sheets appear empty in codon table XLSX: {path}")
        sheet_name = chosen
    else:
        if sheet_name not in sheets:
            raise InputFormatError(
                f"Sheet {sheet_name!r} not found in {path}. Available: {sheets}"
            )

    df = xl.parse(sheet_name)
    if df is None or df.empty:
        raise InputFormatError(f"Codon table sheet {sheet_name!r} is empty in {path}")

    # Normalize column names for detection
    cols = {_clean_colname(c): c for c in df.columns}

    # Detect long format
    aa_col = None
    codon_col = None
    frac_col = None

    for key in ("aa", "aminoacid", "aminoacidletter", "aminoacid1", "amin oacid"):
        if key in cols:
            aa_col = cols[key]
            break
    if aa_col is None:
        # common variants
        for key in cols:
            if key in ("aminoacid", "aminoacid1", "aa"):
                aa_col = cols[key]
                break

    for key in ("codon",):
        if key in cols:
            codon_col = cols[key]
            break

    # fraction/frequency column variants
    for key in ("fraction", "frequency", "usage", "frac"):
        if key in cols:
            frac_col = cols[key]
            break

    # Detect wide format if no codon col
    codons_col = None
    if codon_col is None:
        for key in ("codons", "codonlist"):
            if key in cols:
                codons_col = cols[key]
                break

    table: CodonTable = {}

    if aa_col is not None and codon_col is not None:
        # Long format
        for _, row in df.iterrows():
            if pd.isna(row.get(aa_col)) or pd.isna(row.get(codon_col)):
                continue
            aa = _norm_aa(row[aa_col])
            codon = _norm_codon(row[codon_col])

            if frac_col is None or pd.isna(row.get(frac_col)):
                fr = 1.0
            else:
                try:
                    fr = float(row[frac_col])
                except Exception:
                    raise InputFormatError(
                        f"Invalid fraction value {row[frac_col]!r} for AA {aa}, codon {codon} "
                        f"in {path} sheet {sheet_name!r}"
                    )

            table.setdefault(aa, []).append((codon, fr))

        # Normalize per AA
        for aa in list(table.keys()):
            table[aa] = _normalize_fractions(table[aa])

        if not table:
            raise InputFormatError(
                f"Parsed no codons from long-format table in {path} sheet {sheet_name!r}"
            )

        return table

    if aa_col is not None and codons_col is not None:
        # Wide format
        for _, row in df.iterrows():
            if pd.isna(row.get(aa_col)) or pd.isna(row.get(codons_col)):
                continue
            aa = _norm_aa(row[aa_col])
            codons = _parse_codons_cell(row[codons_col])
            if not codons:
                continue
            # uniform weights unless a fraction column exists (rare in wide format)
            rows = [(c, 1.0) for c in codons]
            table.setdefault(aa, []).extend(rows)

        for aa in list(table.keys()):
            # deduplicate while preserving order
            seen = set()
            dedup = []
            for c, fr in table[aa]:
                if c in seen:
                    continue
                seen.add(c)
                dedup.append((c, fr))
            table[aa] = _normalize_fractions(dedup)

        if not table:
            raise InputFormatError(
                f"Parsed no codons from wide-format table in {path} sheet {sheet_name!r}"
            )
        return table

    # Could not detect format
    raise InputFormatError(
        f"Could not detect codon table format in {path} sheet {sheet_name!r}. "
        f"Expected columns for long format (AA + Codon + Fraction) or wide format (Amino Acid + Codons). "
        f"Found columns: {list(df.columns)}"
    )
