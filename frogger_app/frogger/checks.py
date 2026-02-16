from Bio.Seq import Seq


def _max_homopolymer(seq: str) -> int:
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


def _find_all_overlapping(haystack: str, needle: str):
    hits = []
    start = 0
    while True:
        i = haystack.find(needle, start)
        if i == -1:
            break
        hits.append(i)
        start = i + 1
    return hits


def _context(seq: str, start: int, end: int, flank: int = 20) -> str:
    a = max(0, start - flank)
    b = min(len(seq), end + flank)
    return seq[a:b]


def check_construct(
    construct: dict,
    codonopt_cfg: dict,
    final_forbidden_motifs=None,
    require_no_internal_stops: bool = True,
    require_len_multiple_of_3: bool = True,
):
    """
    Checks are performed on the POST-LIGATION ORF (core fragments only).
    We scan for forbidden motifs derived from codonopt avoid_motifs plus any extra final motifs.
    """
    seq = construct["core_orf_seq"].upper()
    fails = []
    hit_rows = []

    if require_len_multiple_of_3 and (len(seq) % 3 != 0):
        fails.append("len_not_multiple_of_3")

    aa = str(Seq(seq).translate(to_stop=False))
    if require_no_internal_stops and ("*" in aa[:-1]):  # allow terminal stop
        fails.append("internal_stop_codon")

    # homopolymer
    hp = _max_homopolymer(seq)
    max_hp_allowed = int(codonopt_cfg.get("max_homopolymer", 5))
    if hp > max_hp_allowed:
        fails.append(f"homopolymer>{max_hp_allowed}")

    # GC
    gc = (seq.count("G") + seq.count("C")) / max(1, len(seq))
    gc_min = codonopt_cfg.get("gc_min", None)
    gc_max = codonopt_cfg.get("gc_max", None)
    if gc_min is not None and gc < float(gc_min):
        fails.append("gc_below_min")
    if gc_max is not None and gc > float(gc_max):
        fails.append("gc_above_max")

    # forbidden motifs
    forbidden = []
    forbidden += [m.upper() for m in (codonopt_cfg.get("avoid_motifs", []) or []) if str(m).strip()]
    forbidden += [m.upper() for m in (final_forbidden_motifs or []) if str(m).strip()]

    for motif in forbidden:
        for start in _find_all_overlapping(seq, motif):
            end = start + len(motif)
            hit_rows.append(
                {
                    "motif": motif,
                    "start": start,
                    "end": end,
                    "context": _context(seq, start, end, flank=20),
                }
            )

    if hit_rows:
        fails.append("forbidden_motif_or_site")

    return {
        "passes_all": len(fails) == 0,
        "fail_reasons": fails,
        "gc_fraction": gc,
        "max_homopolymer": hp,
        "aa_len": len(aa),
        "hits": hit_rows,
    }
