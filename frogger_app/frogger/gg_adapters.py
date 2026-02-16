from __future__ import annotations

from typing import Dict, List, Optional

from frogger.overhangs import revcomp


ENZYMES = {
    "BSAI": "GGTCTC",
    "BBSI": "GAAGAC",
    "BSMBI": "CGTCTC",
}


def enzyme_site(enzyme: str) -> str:
    key = str(enzyme).upper().replace(" ", "")
    if key not in ENZYMES:
        raise ValueError(f"Unsupported enzyme '{enzyme}'. Supported: {sorted(ENZYMES)}")
    return ENZYMES[key]


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
) -> Dict[str, str]:
    site = enzyme_site(enzyme)
    left = f"{primer_left}{site}{spacer_left}{left_overhang}"
    right_site = revcomp(site)
    right = f"{right_overhang}{spacer_right}{right_site}{primer_right}"
    return {
        "primer_left": primer_left,
        "primer_right": primer_right,
        "left_overhang": left_overhang,
        "right_overhang": right_overhang,
        "left_adapter": left,
        "right_adapter": right,
        "cloning_seq": f"{left}{core_seq}{right}",
    }


def build_cloning_sequences_for_fragments(
    fragments: List[Dict],
    enzyme: str,
    overhangs_by_pos: Dict[int, Dict[str, str]],
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
        )
        out.append({**f, **built})
    return out
