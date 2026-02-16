from __future__ import annotations

from pathlib import Path
from typing import Optional

from frogger.config import load_config
from frogger.pipeline import run_pipeline


def run(config_path: str, outdir: Optional[str] = None) -> None:
    """
    Shuffle workflow: calls the existing FROGGER AA-first pipeline.
    """
    cfg = load_config(config_path)

    if not outdir:
        raise ValueError(
            "No outdir provided. Your current config.yaml does not define outputs.outdir, "
            "so please pass --outdir."
        )

    run_pipeline(cfg, Path(outdir))
