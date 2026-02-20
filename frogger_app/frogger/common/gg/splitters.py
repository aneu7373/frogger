from __future__ import annotations

from typing import Any, Dict, List, Sequence, Tuple


def split_sequence_by_cut_points(
    seq: str,
    cut_points: List[int],
    enforce_frame: bool = True,
) -> List[Tuple[int, int, str]]:
    """
    cut_points: 0-based cut points into seq (nucleotide indices for CDS)
    Example: [300,600] => [0:300],[300:600],[600:end]
    """
    if any(cp < 0 or cp > len(seq) for cp in cut_points):
        raise ValueError("Cut points must be within [0, len(seq)].")
    if sorted(cut_points) != list(cut_points):
        raise ValueError("Cut points must be sorted ascending.")
    if enforce_frame and any((cp % 3) != 0 for cp in cut_points):
        raise ValueError("Cut points must be multiples of 3 when enforce_frame=true.")

    points = [0] + list(cut_points) + [len(seq)]
    frags = []
    for a, b in zip(points[:-1], points[1:]):
        if b <= a:
            raise ValueError("Invalid cut points (non-increasing).")
        frags.append((a, b, seq[a:b]))
    return frags


def aa_cut_points_to_nt_cut_points(cds: str, aa_cut_points: Sequence[int]) -> List[int]:
    """
    Convert AA cut points (e.g. [116, 236]) into nt cut points (0-based).
    Validates against CDS length.
    """
    seq = (cds or "").strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
    if len(seq) % 3 != 0:
        raise ValueError(f"CDS length must be multiple of 3 (len={len(seq)}).")
    aa_len = len(seq) // 3

    cps = [int(x) for x in (aa_cut_points or [])]
    for x in cps:
        if x <= 0 or x >= aa_len:
            raise ValueError(f"Invalid AA cut point {x} for aa_len={aa_len}. Must be within (0, aa_len).")

    cps = sorted(set(cps))
    return [3 * x for x in cps]

def split_aa_by_cut_points(
    aa_seq: str,
    aa_cut_points: List[int],
) -> List[Tuple[int, int, str]]:
    """
    Split an amino-acid sequence by AA cut points (0-based indices).
    Example: aa_cut_points=[116,236] => [0:116], [116:236], [236:end]
    """
    seq = (aa_seq or "").strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
    # For AA splitting, enforce_frame must be False; we just call the primitive.
    return split_sequence_by_cut_points(seq, list(aa_cut_points), enforce_frame=False)

def split_cds_to_fragment_dicts(
    cds: str,
    *,
    aa_cut_points: Sequence[int],
    enforce_frame: bool = True,
) -> List[Dict[str, Any]]:
    """
    Convenience: split CDS by AA cut points and return fragment dicts compatible with gg_adapters.
    """
    nt_cps = aa_cut_points_to_nt_cut_points(cds, aa_cut_points)
    frags_nt = split_sequence_by_cut_points(cds, nt_cps, enforce_frame=enforce_frame)

    out: List[Dict[str, Any]] = []
    for i, (a, b, core) in enumerate(frags_nt, start=1):
        out.append(
            {
                "frag_index": i,
                "core_seq": core,
                "nt_start0": a,
                "nt_end0": b,
                "aa_start0": a // 3,
                "aa_end0": b // 3,
            }
        )
    return out

from typing import Any, Dict, List, Sequence, Tuple

def split_cds_by_aa_cuts_with_4bp_overhang_overlap(
    cds: str,
    *,
    aa_cut_points: Sequence[int],
    window_offset_by_pos: Dict[int, int],
) -> List[Dict[str, Any]]:
    """
    Fix A splitter.

    For each internal junction at AA cut point cp (boundary before aa[cp]),
    we select an overhang window inside the 6 nt of (codon[cp] + codon[cp+1]):

        overhang_start_nt = 3*cp + window_offset
        overhang = cds[overhang_start_nt : overhang_start_nt+4]

    We then generate overlapping fragments so that:
      - left fragment ENDS after the overhang (includes the 4 nt)
      - right fragment STARTS at the overhang (also includes the same 4 nt)

    That makes the Golden Gate ligation *seamless* (no insertion), because the 4 nt
    are an overlap, not an added block.
    """
    seq = (cds or "").strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
    if len(seq) % 3 != 0:
        raise ValueError(f"CDS length must be multiple of 3 (len={len(seq)}).")

    aa_len = len(seq) // 3
    cps = [int(x) for x in (aa_cut_points or [])]
    cps = sorted(cps)

    for x in cps:
        if x <= 0 or x >= aa_len:
            raise ValueError(f"Invalid AA cut point {x} for aa_len={aa_len}. Must be within (0, aa_len).")

    # Compute junction start indices (where the overlap begins)
    # pos_index is 1..len(cps)
    junction_start_nt: List[int] = []
    for pos_index, cp in enumerate(cps, start=1):
        off = int(window_offset_by_pos.get(pos_index, 0))
        if off not in (0, 1, 2):
            raise ValueError(f"window_offset must be 0/1/2, got {off} at pos_index={pos_index}")
        s = 3 * cp + off
        if s < 0 or s + 4 > len(seq):
            raise ValueError(f"Computed overhang start out of range: pos_index={pos_index} start_nt={s}")
        junction_start_nt.append(s)

    # Fragment starts: 0, then each junction_start_nt
    starts = [0] + junction_start_nt
    # Fragment ends: for all but last, end = next_start + 4 (include the overlap)
    ends = []
    for i in range(len(starts) - 1):
        ends.append(starts[i + 1] + 4)
    ends.append(len(seq))

    out: List[Dict[str, Any]] = []
    for i, (a, b) in enumerate(zip(starts, ends), start=1):
        if b <= a:
            raise ValueError("Invalid computed fragment ranges (non-increasing).")
        core = seq[a:b]
        out.append(
            {
                "frag_index": i,
                "core_seq": core,
                "nt_start0": a,
                "nt_end0": b,
                "aa_start0": a // 3,
                "aa_end0": (b // 3),  # note: may not align cleanly at internal overlaps; used mainly for reporting
            }
        )
    return out
