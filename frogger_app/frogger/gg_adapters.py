from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

from frogger.overhangs import revcomp


@dataclass(frozen=True)
class EnzymeDef:
    name: str
    site: str          # recognition site as written 5'->3' on the "left" side
    top_cut: int       # a in NEB notation: site(a/b)
    bottom_cut: int    # b in NEB notation: site(a/b)

    @property
    def overhang_len(self) -> int:
        return self.bottom_cut - self.top_cut


# Cut-site definitions from vendor conventions:
# BsaI:  GGTCTC(1/5)
# BbsI:  GAAGAC(2/6)
# BsmBI: CGTCTC(1/5)  (Esp3I is an isoschizomer commonly used interchangeably)
# PaqCI: CACCTGC(4/8)
#
# These yield 4-nt overhangs in each case.
ENZYMES: Dict[str, EnzymeDef] = {
    "BSAI": EnzymeDef("BsaI", "GGTCTC", 1, 5),
    "BBSI": EnzymeDef("BbsI", "GAAGAC", 2, 6),
    "BSMBI": EnzymeDef("BsmBI", "CGTCTC", 1, 5),
    "ESP3I": EnzymeDef("Esp3I", "CGTCTC", 1, 5),
    "PAQCI": EnzymeDef("PaqCI", "CACCTGC", 4, 8),
}


def _get_enzyme(enzyme: str) -> EnzymeDef:
    key = str(enzyme).upper().replace(" ", "")
    if key not in ENZYMES:
        raise ValueError(f"Unsupported enzyme '{enzyme}'. Supported: {sorted(ENZYMES)}")
    return ENZYMES[key]


def _offset_filler(a: int) -> str:
    """
    The offset filler bases between the recognition site and the start of the 4-nt overhang.
    These bases do NOT matter for overhang identity, but MUST exist to align with the cut geometry.

    We choose 'A'*a deterministically.
    """
    return "A" * int(a)


def add_adapters(fragments: List[Dict], gg_cfg: List[Dict]) -> List[Dict]:
    """
    Legacy adapter concatenation (kept for backward compatibility).
    """
    out = []
    for f in fragments:
        idx = int(f["frag_index"]) - 1
        if idx < 0 or idx >= len(gg_cfg):
            raise ValueError(f"No Golden Gate adapter specified for frag_index={f['frag_index']}")
        left = gg_cfg[idx]["left_adapter"]
        right = gg_cfg[idx]["right_adapter"]
        out.append({**f, "cloning_seq": f"{left}{f['core_seq']}{right}"})
    return out


def build_cloning_seq(
    core_seq: str,
    enzyme: str,
    left_overhang: str,
    right_overhang: str,
    primer_left: str = "",
    primer_right: str = "",
    spacer_left: str = "",
    spacer_right: str = "",
    embed_left_overhang_in_core: bool = False,
    embed_right_overhang_in_core: bool = False,
) -> Dict[str, str]:
    """
    Build a Golden Gate cloning block.

    IMPORTANT: This implementation follows the user's "template-strand sequential ligation" convention.

    - left_overhang is the desired 4-nt sticky end sequence written 5'->3' on the TOP (template) strand.
      (This is the sequence that appears at the start of the assembled contig at that junction.)

    - right_overhang is the NEXT fragment's left_overhang label (i.e., the junction label on the TOP strand).
      For sequential ligation, the RIGHT end of the current fragment must expose the COMPLEMENT of that label
      on the opposite strand. Example:
          next left = TGAT  ==> current right must expose ACTA (complement) on the bottom strand.

    Type IIS geometry:
    - We insert ez.top_cut filler bases between the recognition site and the overhang region.

    Spacer semantics (safety fix):
    - spacer_left/spacer_right are placed OUTSIDE the recognition site, so they cannot accidentally change
      sticky-end identity if non-empty.

    Right-end construction (core fix for your convention):
    - To make the *physical* 5' overhang on the bottom strand equal complement(right_overhang),
      we must write reverse(right_overhang) on the TOP strand immediately upstream of the right site_rc.
      (Not revcomp.)
    """
    ez = _get_enzyme(enzyme)
    if ez.overhang_len != 4:
        raise ValueError(f"{ez.name} configured with overhang_len={ez.overhang_len}, expected 4.")

    left_overhang = str(left_overhang).upper().strip()
    right_overhang = str(right_overhang).upper().strip()

    if len(left_overhang) != 4 or len(right_overhang) != 4:
        raise ValueError("Overhangs must be exactly 4 bp each.")
    if any(b not in "ACGT" for b in left_overhang + right_overhang):
        raise ValueError("Overhangs must contain only A/C/G/T.")

    site = ez.site
    site_rc = revcomp(site)

    # Enzyme offset filler (length = top_cut)
    filler = _offset_filler(ez.top_cut)
    filler_rc = revcomp(filler)

    # LEFT end (top strand):
    # If embed_left_overhang_in_core=True, the core already begins with the 4 nt that define the sticky end.
    if embed_left_overhang_in_core:
        left_adapter = f"{primer_left}{site}{spacer_left}{filler}"
    else:
        left_adapter = f"{primer_left}{site}{spacer_left}{filler}{left_overhang}"

    # RIGHT end:
    # If embed_right_overhang_in_core=True, the core already ENDS with revcomp(right_overhang) on the top strand.
    right_overhang_rc = revcomp(right_overhang)
    if embed_right_overhang_in_core:
        right_adapter = f"{filler_rc}{spacer_right}{site_rc}{primer_right}"
    else:
        right_adapter = f"{right_overhang_rc}{filler_rc}{spacer_right}{site_rc}{primer_right}"

    cloning_seq = f"{left_adapter}{core_seq}{right_adapter}"

    return {
        "primer_left": primer_left,
        "primer_right": primer_right,
        "left_overhang": left_overhang,
        "right_overhang": right_overhang,
        "left_offset_filler": filler,
        "right_offset_filler_rc": filler_rc,
        "right_overhang_rc": right_overhang_rc,
        "left_adapter": left_adapter,
        "right_adapter": right_adapter,
        "cloning_seq": cloning_seq,
        "enzyme_site": site,
        "enzyme_site_rc": site_rc,
    }


def build_cloning_sequences_for_fragments(
    fragments: List[Dict],
    enzyme: str,
    overhangs_by_pos: Dict[int, Dict[str, str]],
    n_fragments: int,
    primers_by_pos: Optional[Dict[int, Dict[str, str]]] = None,
    spacer_left: str = "",
    spacer_right: str = "",
) -> List[Dict]:
    primers_by_pos = primers_by_pos or {}
    out = []
    for f in fragments:
        idx = int(f["frag_index"])
        if idx not in overhangs_by_pos:
            raise ValueError(f"Missing overhang assignment for frag_index={idx}")
        oh = overhangs_by_pos[idx]
        pr = primers_by_pos.get(idx, {})
        built = build_cloning_seq(
            core_seq=f["core_seq"],
            enzyme=enzyme,
            left_overhang=oh["left"],
            right_overhang=oh["right"],
            primer_left=pr.get("left", ""),
            primer_right=pr.get("right", ""),
            spacer_left=spacer_left,
            spacer_right=spacer_right,
            embed_left_overhang_in_core=(idx > 1),
            embed_right_overhang_in_core=(idx < n_fragments),
        )
        out.append({**f, **built})
    return out

def _comp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))


def _extract_left_label_from_cloning_seq(cloning_seq: str, enzyme: str) -> str:
    ez = _get_enzyme(enzyme)
    site = ez.site
    i = cloning_seq.find(site)
    if i < 0:
        raise ValueError(f"Left site {site} not found in cloning_seq")
    start = i + len(site) + ez.top_cut
    return cloning_seq[start:start + 4].upper()


def _extract_right_top_written_overhang_block(cloning_seq: str, enzyme: str, spacer_right: str = "") -> str:
    """
    Extract the 4-bp OVERHANG BLOCK that was written on the TOP strand immediately upstream
    of the (filler_rc + spacer_right + site_rc) region.

    IMPORTANT:
      right_adapter is built as:
        [overhang_block (4)] + [filler_rc (top_cut)] + [spacer_right] + [site_rc] + ...
      Therefore the 4 bases immediately before site_rc are NOT the overhang block when top_cut>0.
    """
    ez = _get_enzyme(enzyme)
    site_rc = revcomp(ez.site)
    j = cloning_seq.rfind(site_rc)
    if j < 0:
        raise ValueError(f"Right site_rc {site_rc} not found in cloning_seq")

    a = int(ez.top_cut)
    s = len(spacer_right or "")
    end = j - (a + s)          # end of the 4bp block
    start = end - 4
    if start < 0:
        raise ValueError("Cloning seq too short to extract right overhang block")
    return cloning_seq[start:end].upper()


def validate_physical_golden_gate_chain(
    frags: List[Dict[str, str]],
    enzyme: str,
) -> None:
    """
    Validate adjacency under the CANONICAL label convention.

    Definitions:
      - left_label(next) is the canonical 5'->3' sticky-end label for the LEFT end of the next fragment.
      - physical_right_bottom_5p(current) is the 5'->3' sequence of the protruding sticky end exposed
        on the BOTTOM strand at the RIGHT end of the current fragment after Type IIS digestion.

    Canonical rule:
        physical_right_bottom_5p(current) == left_label(next)

    Note:
      This function assumes each fragment dict contains:
        - "cloning_seq"
        - optionally "spacer_right" (or we assume "")
      and that `frags` is already in assembly order.
    """
    if frags is None or len(frags) < 2:
        return

    for k in range(len(frags) - 1):
        a = frags[k]
        b = frags[k + 1]

        if "cloning_seq" not in a or "cloning_seq" not in b:
            raise ValueError("validate_physical_golden_gate_chain: missing cloning_seq in fragment dict")

        seq_a = a["cloning_seq"]
        seq_b = b["cloning_seq"]

        # Canonical left label from the next fragment
        left_b = _extract_left_label_from_cloning_seq(seq_b, enzyme)

        # Extract the 4bp block written upstream of (filler_rc + spacer_right + site_rc)
        top_written_overhang_block_a = _extract_right_top_written_overhang_block(
            seq_a,
            enzyme,
            spacer_right=a.get("spacer_right", ""),
        )

        # Physical 5'->3' sticky end on the protruding bottom strand
        physical_right_bottom_a = revcomp(top_written_overhang_block_a)

        expected = left_b

        if physical_right_bottom_a != expected:
            raise ValueError(
                "Golden Gate physical adjacency FAILED under CANONICAL label convention:\n"
                f"  fragment[{k}] -> fragment[{k+1}]\n"
                f"  next_left_label(5'->3')          = {left_b}\n"
                f"  expected_right_bottom_5p         = {expected}\n"
                f"  current_top_written_overhangblk  = {top_written_overhang_block_a}\n"
                f"  current_right_bottom_5p(observed)= {physical_right_bottom_a}\n"
            )
