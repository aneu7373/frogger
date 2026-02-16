import subprocess
from pathlib import Path
from typing import Optional, Tuple, List


def effective_seed(user_seed: Optional[int], default_seed: int) -> int:
    return int(default_seed if user_seed is None else user_seed)


def _join_bar(items) -> str:
    items = [str(x).strip() for x in (items or []) if str(x).strip()]
    return "|".join(items)


def _build_codonopt_cmd(c: dict, outdir: Path, seed: int, sequence_path: str) -> List[str]:
    cmd = ["python", "-m", "codonopt.main"]

    # codonopt supports --sequence for FASTA input; we pass fragment FASTA here
    cmd += ["--sequence", sequence_path]

    if not c.get("codon_table"):
        raise ValueError("codonopt.codon_table is required.")
    cmd += ["--codon-table", c["codon_table"]]
    if c.get("codon_table_sheet"):
        cmd += ["--codon-table-sheet", c["codon_table_sheet"]]

    cmd += ["--out", str(outdir)]
    if c.get("log_file"):
        cmd += ["--log-file", c["log_file"]]
    if bool(c.get("verbose", False)):
        cmd += ["--verbose"]

    cmd += ["--avoid-codons", _join_bar(c.get("avoid_codons"))]
    cmd += ["--avoid-motifs", _join_bar(c.get("avoid_motifs"))]

    if c.get("gc_min") is not None:
        cmd += ["--gc-min", str(float(c["gc_min"]))]

    if c.get("gc_max") is not None:
        cmd += ["--gc-max", str(float(c["gc_max"]))]

    cmd += ["--max-homopolymer", str(int(c.get("max_homopolymer", 5)))]
    cmd += ["--min-codon-fraction", str(float(c.get("min_codon_fraction", 0.05)))]
    cmd += ["--n", str(int(c.get("n", 1)))]
    cmd += ["--seed", str(int(seed))]

    opt = c.get("optimization", {}) or {}
    cmd += ["--optimization-mode", str(opt.get("mode", "kleinbub"))]
    cmd += ["--max-tries-per-replicate", str(int(opt.get("max_tries_per_replicate", 25)))]
    cmd += ["--kleinbub-search-limit", str(int(opt.get("kleinbub_search_limit", 200000)))]
    cmd += ["--backtrack-window", str(int(opt.get("backtrack_window", 10)))]

    ck = c.get("cryptkeeper", {}) or {}
    if bool(ck.get("enable", False)):
        cmd += ["--cryptkeeper-enable"]
        cmd += ["--cryptkeeper-rbs-score-cutoff", str(float(ck.get("rbs_score_cutoff", 500.0)))]
        cmd += ["--cryptkeeper-threads", str(int(ck.get("threads", 1)))]
        cmd += ["--cryptkeeper-fail-score", str(float(ck.get("fail_score", 6000.0)))]
        cmd += ["--cryptkeeper-ignore-first-nt", str(int(ck.get("ignore_first_nt", 30)))]
        cmd += ["--cryptkeeper-ignore-last-nt", str(int(ck.get("ignore_last_nt", 30)))]
        cmd += ["--cryptkeeper-pool-factor", str(float(ck.get("pool_factor", 3.0)))]
        cmd += ["--cryptkeeper-max-pool", str(int(ck.get("max_pool", 200)))]

    return cmd


def run_codonopt_for_fragment(
    codonopt_cfg: dict,
    outdir: Path,
    seed: int,
    fragment_id: str,
    aa_seq: str,
) -> Tuple[str, Path, Path]:
    """
    Run codonopt on a single AA fragment (written as a temporary FASTA).
    Returns (optimized_cds, fasta_path, metrics_path).
    """
    frag_out = outdir / "codonopt_fragments" / fragment_id
    frag_out.mkdir(parents=True, exist_ok=True)

    tmp_fa = frag_out / "input_fragment.faa"
    with open(tmp_fa, "w") as f:
        f.write(f">{fragment_id}\n{aa_seq}\n")

    cmd = _build_codonopt_cmd(codonopt_cfg, frag_out, seed, str(tmp_fa))

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        msg = (
            f"codonopt failed for fragment '{fragment_id}' (exit {proc.returncode}).\n"
            f"STDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}\n"
        )
        raise RuntimeError(msg)

    fasta = frag_out / "optimized_sequences.fasta"
    metrics = frag_out / "metrics.tsv"
    if not fasta.exists():
        raise FileNotFoundError(f"Expected codonopt output missing: {fasta}")
    if not metrics.exists():
        raise FileNotFoundError(f"Expected codonopt output missing: {metrics}")

    optimized = None
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                continue
            optimized = (optimized or "") + line.strip()
    if not optimized:
        raise RuntimeError(f"No optimized sequence found in {fasta}")

    return optimized, fasta, metrics
