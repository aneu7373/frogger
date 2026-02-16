# frogger/common/qc.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple


_DNA = set("ACGT")


def _norm_dna(s: str) -> str:
    s = (s or "").upper().replace("U", "T")
    return "".join([c for c in s if c in _DNA])


def revcomp(seq: str) -> str:
    seq = _norm_dna(seq)
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(comp[b] for b in reversed(seq))


def max_homopolymer_run(seq: str) -> Tuple[int, str]:
    """
    Returns (max_run_length, base) e.g. (7, 'A').
    """
    s = _norm_dna(seq)
    if not s:
        return 0, ""
    best_len = 1
    best_base = s[0]
    cur_len = 1
    cur_base = s[0]
    for b in s[1:]:
        if b == cur_base:
            cur_len += 1
        else:
            if cur_len > best_len:
                best_len, best_base = cur_len, cur_base
            cur_base = b
            cur_len = 1
    if cur_len > best_len:
        best_len, best_base = cur_len, cur_base
    return best_len, best_base


@dataclass(frozen=True)
class Site:
    name: str
    site: str  # recognition sequence (DNA)


def _find_all(haystack: str, needle: str) -> List[int]:
    """All start indices of needle in haystack (non-overlapping is NOT enforced)."""
    if not needle:
        return []
    out = []
    i = 0
    while True:
        j = haystack.find(needle, i)
        if j < 0:
            break
        out.append(j)
        i = j + 1
    return out


def qc_scan_nt(seq: str, qc_cfg: dict) -> Dict[str, object]:
    """
    Minimal QC:
      - max_homopolymer
      - forbidden_motifs (and their revcomps)
      - internal_sites (and their revcomps)
    Returns:
      { "pass": bool, "failures": [str...], "metrics": {...} }
    """
    qc_cfg = qc_cfg or {}
    s = _norm_dna(seq)

    failures: List[str] = []
    metrics: Dict[str, object] = {}

    # homopolymers
    mh = qc_cfg.get("max_homopolymer")
    if mh is not None:
        mh = int(mh)
        run_len, base = max_homopolymer_run(s)
        metrics["max_homopolymer"] = run_len
        metrics["max_homopolymer_base"] = base
        if run_len > mh:
            failures.append(f"HOMOPOLYMER>{base}:{run_len}")

    # motifs
    motifs = qc_cfg.get("forbidden_motifs", []) or []
    motif_hits = 0
    for m in motifs:
        m = _norm_dna(m)
        if not m:
            continue
        rc = revcomp(m)
        hits_f = _find_all(s, m)
        hits_r = _find_all(s, rc) if rc != m else []
        if hits_f:
            motif_hits += len(hits_f)
            failures.append(f"MOTIF:{m}:{len(hits_f)}")
        if hits_r:
            motif_hits += len(hits_r)
            failures.append(f"MOTIF_REVCOMP:{m}:{len(hits_r)}")
    metrics["forbidden_motif_hits"] = motif_hits

    # internal restriction sites
    sites_cfg = qc_cfg.get("internal_sites", []) or []
    site_hits = 0
    for entry in sites_cfg:
        # accept dicts {name, site} or tuples
        if isinstance(entry, dict):
            name = str(entry.get("name", "SITE"))
            site = _norm_dna(entry.get("site", ""))
        else:
            name, site = entry  # type: ignore[misc]
            name, site = str(name), _norm_dna(str(site))
        if not site:
            continue
        rc = revcomp(site)
        hits_f = _find_all(s, site)
        hits_r = _find_all(s, rc) if rc != site else []
        if hits_f:
            site_hits += len(hits_f)
            failures.append(f"SITE:{name}:{site}:{len(hits_f)}")
        if hits_r:
            site_hits += len(hits_r)
            failures.append(f"SITE_REVCOMP:{name}:{site}:{len(hits_r)}")
    metrics["internal_site_hits"] = site_hits

    return {"pass": (len(failures) == 0), "failures": failures, "metrics": metrics}
