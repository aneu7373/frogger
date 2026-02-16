from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

from frogger import __version__


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="frogger",
        description="FROGGER: Golden Gate–ready fragment library generator (shuffle + mutate).",
    )
    p.add_argument("--version", action="version", version=f"frogger {__version__}")

    sub = p.add_subparsers(dest="command", required=True)

    shuf = sub.add_parser("shuffle", help="Homolog shuffling / fragment library generation.")
    shuf.add_argument("--config", required=True, help="Path to YAML config file.")
    shuf.add_argument("--outdir", required=True, help="Output directory (required).")

    mut = sub.add_parser("mutate", help="Mutagenesis library generation (SLOGGER port target).")
    mut.add_argument("--config", required=True, help="Path to YAML config file.")
    mut.add_argument("--outdir", required=True, help="Output directory (required).")

    run = sub.add_parser("run", help="(Legacy) Alias for `shuffle` for backward compatibility.")
    run.add_argument("--config", required=True, help="Path to YAML config file.")
    run.add_argument("--outdir", required=True, help="Output directory (required).")
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = _build_parser().parse_args(argv)

    cfg_path = Path(args.config)
    if not cfg_path.exists():
        print(f"ERROR: config file not found: {cfg_path}", file=sys.stderr)
        return 2

    if args.command in ("shuffle", "run"):
        from frogger.modules.shuffle.pipeline import run as run_shuffle
        run_shuffle(config_path=args.config, outdir=args.outdir)
        return 0

    if args.command == "mutate":
        from frogger.modules.mutate.pipeline import run as run_mutate
        run_mutate(config_path=args.config, outdir=args.outdir)
        return 0

    print(f"ERROR: Unknown command: {args.command}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
