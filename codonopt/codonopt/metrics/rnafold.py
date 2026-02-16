
# codonopt/metrics/rnafold.py
"""
RNAfold utilities (disabled by default).

This module provides optional helpers to compute RNA minimum free energy (MFE)
using the ViennaRNA RNAfold binary. The main program should *not* require
RNAfold to run. If RNAfold is disabled (default) or unavailable, the helpers
either return None or raise a clear RuntimeError when explicitly requested.

Usage pattern in the rest of the code:
    - Pass enable_rnafold=False (default) to skip all RNAfold calls.
    - If enable_rnafold=True, we try to run RNAfold and handle missing binary.

CLI should expose a flag like: --enable-rnafold (default: off)
"""

from __future__ import annotations

import os
import shutil
import subprocess
from typing import Optional, Tuple


def _rnafold_available() -> bool:
    """Return True if RNAfold binary is on PATH."""
    return shutil.which("RNAfold") is not None


def _run_rnafold(seq: str, extra_args: Optional[list[str]] = None, timeout: int = 30) -> Tuple[str, float]:
    """
    Run RNAfold on a single RNA sequence string and return (structure, mfe_kcal_per_mol).

    Parameters
    ----------
    seq : str
        RNA sequence (A/C/G/U). If your pipeline gives DNA, convert T->U first.
    extra_args : list[str] | None
        Additional CLI arguments to pass to RNAfold, e.g., ["--noPS"].
    timeout : int
        Seconds to wait before killing the process.

    Returns
    -------
    (structure, mfe)

    Raises
    ------
    RuntimeError if RNAfold is missing or returns a non-zero exit code.
    """
    if not _rnafold_available():
        raise RuntimeError("RNAfold binary not found in PATH. Install ViennaRNA or disable RNAfold.")

    args = ["RNAfold", "--noPS"]
    if extra_args:
        args.extend(extra_args)

    # RNAfold reads the sequence from stdin and prints two lines:
    # 1) the sequence
    # 2) the structure + MFE in the form: "((((....)))) (-12.34)"
    try:
        proc = subprocess.run(
            args,
            input=(seq + "\n").encode(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=timeout,
        )
    except FileNotFoundError:
        raise RuntimeError("RNAfold binary not found in PATH. Install ViennaRNA or disable RNAfold.")
    except subprocess.TimeoutExpired:
        raise RuntimeError("RNAfold timed out.")

    if proc.returncode != 0:
        raise RuntimeError(f"RNAfold failed (exit {proc.returncode}): {proc.stderr.decode(errors='ignore').strip()}")

    out = proc.stdout.decode(errors="ignore").strip().splitlines()
    if len(out) < 2:
        raise RuntimeError(f"Unexpected RNAfold output: {out!r}")

    # Example second line: "....((....)).... (-7.80)"
    line2 = out[1]
    try:
        struct, paren = line2.rsplit(" ", 1)
        mfe = float(paren.strip("() "))
    except Exception:
        # Fallback: try to find a float anywhere in the line
        import re

        m = re.search(r"([-+]?\d+(\.\d+)?)", line2)
        if not m:
            raise RuntimeError(f"Could not parse MFE from RNAfold output line: {line2!r}")
        struct = line2.split()[0]
        mfe = float(m.group(1))

    return struct.strip(), mfe


def rnafold_mfe(rna_seq: str, enable: bool = False) -> Optional[float]:
    """
    Return the MFE (kcal/mol) for the full RNA sequence if enabled, else None.

    - If enable=False (default), returns None.
    - If enable=True and RNAfold is present, returns a float.
    - If enable=True and RNAfold is missing/fails, raises RuntimeError.
    """
    if not enable:
        return None
    _, mfe = _run_rnafold(rna_seq.replace("T", "U").upper())
    return mfe


def rnafold_mfe_5p(rna_seq: str, window: int = 50, enable: bool = False) -> Optional[float]:
    """
    Return the MFE (kcal/mol) for the 5' window (default 50 nt) if enabled, else None.

    - If enable=False (default), returns None.
    - If enable=True and RNAfold is present, returns a float.
    - If enable=True and RNAfold is missing/fails, raises RuntimeError.
    """
    if not enable:
        return None
    subseq = rna_seq[:window].replace("T", "U").upper()
    _, mfe = _run_rnafold(subseq)
    return mfe
