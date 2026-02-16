# frogger/common/repair.py
from __future__ import annotations

from typing import Callable, Dict, List, Tuple

from frogger.common.qc import qc_scan_nt


def repair_fragments_loop(
    fragments_cds: List[Dict],
    *,
    qc_cfg: dict,
    enable: bool,
    max_rounds: int,
    # You inject these from the pipeline so we don't depend on internal signatures here:
    codonopt_fragment_runner: Callable[[Dict], Tuple[str, str]],
    # returns (new_nt_core_seq, codonopt_version)
    aa_translate: Callable[[str], str],
    patch_forced_codons: Callable[[str, int, Dict[int, str]], str] | None = None,
    # if provided, will be applied after codonopt
) -> Tuple[List[Dict], List[Dict]]:
    """
    fragments_cds: list of dicts with at least:
      - fragment_id
      - aa_seq
      - aa_len
      - forced_codons_local (optional dict[int->codon])
      - core_seq (nt)
    Returns:
      (updated_fragments_cds, report_rows)
    """
    if not enable:
        return fragments_cds, []

    max_rounds = max(1, int(max_rounds))
    qc_cfg = qc_cfg or {}

    # Cache original AA to validate preservation
    orig_aa_by_id: Dict[str, str] = {}
    for f in fragments_cds:
        fid = str(f.get("fragment_id") or f.get("id") or "")
        if not fid:
            raise ValueError("repair_fragments_loop: fragment missing fragment_id/id")
        orig_aa_by_id[fid] = str(f.get("aa_seq") or "")

    report_rows: List[Dict] = []

    for rnd in range(1, max_rounds + 1):
        failing: List[Dict] = []
        scan_before: Dict[str, Dict] = {}

        # Scan
        for f in fragments_cds:
            fid = str(f.get("fragment_id"))
            res = qc_scan_nt(str(f.get("core_seq") or ""), qc_cfg)
            scan_before[fid] = res
            if not bool(res["pass"]):
                failing.append(f)

        if not failing:
            break

        # Repair only failing
        for f in failing:
            fid = str(f.get("fragment_id"))
            before = scan_before[fid]
            before_fail = ";".join(before["failures"])

            old_nt = str(f.get("core_seq") or "")
            old_aa = orig_aa_by_id[fid]

            new_nt, codonopt_ver = codonopt_fragment_runner(f)

            # Optionally re-apply forced codons (junction realizability)
            forced_local = f.get("forced_codons_local") or {}
            if patch_forced_codons is not None and forced_local:
                aa_len = int(f.get("aa_len") or len(old_aa))
                new_nt = patch_forced_codons(new_nt, aa_len, forced_local)

            # Validate AA preservation
            new_aa = aa_translate(new_nt)
            aa_ok = (new_aa == old_aa)

            # If AA changed, do NOT replace (keep old) but report
            replaced = False
            if aa_ok and new_nt and new_nt != old_nt:
                f["core_seq"] = new_nt
                replaced = True

            after = qc_scan_nt(str(f.get("core_seq") or ""), qc_cfg)
            after_fail = ";".join(after["failures"])
            status_after = "PASS" if after["pass"] else "FAIL"

            report_rows.append(
                {
                    "fragment_id": fid,
                    "round": rnd,
                    "status_before": "PASS" if before["pass"] else "FAIL",
                    "failures_before": before_fail,
                    "status_after": status_after,
                    "failures_after": after_fail,
                    "aa_ok": int(aa_ok),
                    "nt_changed": int(old_nt != str(f.get("core_seq") or "")),
                    "replaced": int(replaced),
                    "codonopt_version": codonopt_ver or "",
                }
            )

    return fragments_cds, report_rows
