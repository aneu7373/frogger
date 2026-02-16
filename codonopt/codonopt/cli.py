# codonopt/cli.py

import argparse


def get_cli_args():
    parser = argparse.ArgumentParser(
        prog="codonopt",
        description="Codon optimization / back-translation with batch support.",
    )

    # Inputs
    parser.add_argument("--sequence", help="Input FASTA (DNA CDS or protein).")
    parser.add_argument("--batch-table", help="Batch CSV/TSV with per-row parameters.")
    parser.add_argument(
        "--codon-table",
        help="Path to codon table XLSX (used in single-sequence mode).",
    )
    parser.add_argument(
        "--codon-table-sheet",
        default=None,
        help="Optional sheet name in the codon table XLSX (single-sequence mode).",
    )

    # Outputs
    parser.add_argument("--out", default="results", help="Output directory.")
    parser.add_argument("--log-file", default=None, help="Optional log file path.")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")

    # Constraints / knobs
    parser.add_argument("--avoid-codons", default="", help="Codons to avoid, '|' separated.")
    parser.add_argument("--avoid-motifs", default="", help="Motifs to avoid, '|' separated.")
    parser.add_argument("--gc-min", type=float, default=None, help="Minimum GC fraction (0..1).")
    parser.add_argument("--gc-max", type=float, default=None, help="Maximum GC fraction (0..1).")
    parser.add_argument("--max-homopolymer", type=int, default=5, help="Max homopolymer run length.")
    parser.add_argument("--seed", type=int, default=None, help="Seed for reproducibility.")
    parser.add_argument("--n", type=int, default=1, help="Number of sequences to output per input record.")

    # Optimization behavior
    parser.add_argument(
        "--optimization-mode",
        default="kleinbub",
        choices=["kleinbub", "strict"],
        help="Optimization mode. Default kleinbub (bounded backtracking).",
    )
    parser.add_argument(
        "--max-tries-per-replicate",
        type=int,
        default=25,
        help="How many independent restarts to attempt to get ONE passing sequence.",
    )
    parser.add_argument(
        "--kleinbub-search-limit",
        type=int,
        default=200000,
        help="Total search steps per attempt in kleinbub mode (higher = more likely success, slower).",
    )
    parser.add_argument(
        "--backtrack-window",
        type=int,
        default=10,
        help="How far kleinbub can rewind when stuck (higher = more rescue, slower).",
    )

    # v1.1: codon rarity cutoff
    parser.add_argument(
        "--min-codon-fraction",
        type=float,
        default=0.05,
        help="Exclude codons with fraction below this cutoff (default 0.05).",
    )

    # v1.1: CryptKeeper screening
    parser.add_argument(
        "--cryptkeeper-enable",
        action="store_true",
        help="Enable CryptKeeper screening to discard sequences with strong internal translation start sites.",
    )

    # UPDATED DEFAULT: scores in practice are often in the hundreds to thousands
    parser.add_argument(
        "--cryptkeeper-rbs-score-cutoff",
        type=float,
        default=500.0,
        help="CryptKeeper reporting cutoff; sites below this score are ignored (default 500).",
    )

    parser.add_argument(
        "--cryptkeeper-threads",
        type=int,
        default=1,
        help="Threads for CryptKeeper (default 1).",
    )

    # UPDATED DEFAULT: calibrated to typical observed max scores in your runs
    parser.add_argument(
        "--cryptkeeper-fail-score",
        type=float,
        default=3000.0,
        help="Fail a candidate if any internal site has score >= this threshold (default 3000).",
    )

    parser.add_argument(
        "--cryptkeeper-ignore-first-nt",
        type=int,
        default=30,
        help="Ignore predicted start sites starting within the first N nt (default 30).",
    )

    parser.add_argument(
        "--cryptkeeper-ignore-last-nt",
        type=int,
        default=30,
        help="Ignore predicted start sites starting within the last N nt (default 30).",
    )

    # v1.1: pooling strategy (two-phase) for speed
    parser.add_argument(
        "--cryptkeeper-pool-factor",
        type=float,
        default=3.0,
        help="When CryptKeeper is enabled, generate ~n*pool_factor candidates then screen (default 3.0).",
    )

    parser.add_argument(
        "--cryptkeeper-max-pool",
        type=int,
        default=200,
        help="Hard cap on pool size per input record when CryptKeeper is enabled (default 200).",
    )

    return parser.parse_args()
