# COTTAGE: Codon Optimization Tool That Aids Gene Expression

**COTTAGE** is a batch-capable codon optimization / back-translation tool designed for generating many candidate coding sequences under practical constraints (GC bounds, forbidden motifs, forbidden codons, homopolymer limits, codon-usage cutoffs, etc.).

It supports:
- DNA CDS FASTA inputs and protein FASTA inputs
- Batch mode via CSV/TSV with per-row parameters
- Codon tables from Excel (.xlsx) files
- High-yield optimization via a bounded-backtracking mode (“kleinbub”, default)
- Optional filtering of candidates for cryptic internal translation initiation sites using CryptKeeper (v1.2)

Outputs:
- One FASTA containing all passing sequences
- One TSV containing sequence metrics and debug fields

---

## Installation

### Docker (recommended)

    git clone https://github.com/aneu7373/codonopt.git
    cd codonopt
    docker build -t codonopt .

---

## Running

### Single input (FASTA may be DNA CDS or protein)

    docker run --rm \
      -v "$(pwd)/data:/data" \
      codonopt \
      --sequence /data/input.fasta \
      --codon-table /data/Ecoli.xlsx \
      --out /data/results

### Batch mode (recommended)

Batch mode takes a CSV or TSV where each row defines a job. The `sequence` field points to a FASTA file; that FASTA may contain multiple records, and each record will be processed.

Minimal CSV:

    sequence,codon_table
    /data/input.fasta,/data/Ecoli.xlsx

Run:

    docker run --rm \
      -v "$(pwd)/data:/data" \
      codonopt \
      --batch-table /data/batch.csv \
      --out /data/results \
      --verbose

---

## Inputs

### DNA vs protein FASTA

- If the FASTA record contains letters outside A/C/G/T (e.g., amino acids), it is treated as protein.
- Otherwise it is treated as DNA CDS and translated to protein.
- DNA CDS records must have length divisible by 3.

---

## Key constraints

- GC bounds: `gc_min`, `gc_max`
- Avoid motifs: `avoid_motifs` (use `|` separator)
- Avoid codons: `avoid_codons` (use `|` separator)
- Max homopolymer length (default 5): `max_homopolymer`
- Minimum codon fraction (default 0.05): `min_codon_fraction`  
  Any codons for an amino acid with usage fraction < this cutoff are excluded from use.

---

## Optimization modes

- `optimization_mode = kleinbub` (default): bounded backtracking for higher yield under tight constraints
- `optimization_mode = strict`: simpler/faster approach (can fail more often under tight constraints)

Yield tuning knobs:

- `max_tries_per_replicate`: restarts per requested output
- `kleinbub_search_limit`: how hard to search in kleinbub mode
- `backtrack_window`: how far the search can rewind when stuck

---

## Outputs

### FASTA

`optimized_sequences.fasta` contains only passing sequences.

### TSV

`metrics.tsv` contains:

- one row per requested replicate
- sequence metrics (GC, max homopolymer length, etc.)
- debug details (attempt counts, CryptKeeper stats)
- if a replicate fails, `failure_reason` explains why

---

## Codon tables (.xlsx)

Codon tables are loaded from Excel (.xlsx). Two formats are supported.

### Long format (recommended)

| AA | Codon | Fraction |
|----|-------|----------|
| A  | GCT   | 0.18     |
| A  | GCC   | 0.27     |
| ...| ...   | ...      |

Column names can vary (AA/Amino Acid, Codon, Fraction/Frequency/Usage). Fractions are normalized per amino acid if needed.

### Wide format

| Amino Acid | Codons |
|------------|--------|
| A          | GCT,GCC,GCA,GCG |
| C          | TGT,TGC |
| *          | TAA,TAG,TGA |

Wide format assumes uniform weighting among codons.

Multiple sheets are supported:

- batch column: `codon_table_sheet`
- single mode flag: `--codon-table-sheet`

---

# CryptKeeper (v1.1): filtering internal translation initiation sites

## What it does

When enabled, codonopt generates candidates and then uses CryptKeeper to detect potential internal translation initiation sites (TSS/TIS-like signals) in the DNA sequence. Candidates with very strong internal sites can be discarded.

## Important note on score scale (calibration)

In real runs, CryptKeeper scores are often not small numbers (not 0–10). In observed runs, typical max internal scores were in the hundreds to thousands (e.g., ~900 to ~8000). Because of this, a fail threshold like 6.0 will reject everything.

## Defaults (calibrated)

codonopt uses these calibrated defaults:

- `cryptkeeper_rbs_score_cutoff = 500`  
  Only sites with score >= 500 are considered/reportable.
- `cryptkeeper_fail_score = 3000`  
  Reject a candidate only if it contains an internal site with score >= 3000.
- `cryptkeeper_ignore_first_nt = 30`
- `cryptkeeper_ignore_last_nt = 30`

This is intentionally permissive and aims to filter only the most extreme internal sites.

## How pass/fail works

A candidate fails CryptKeeper only if:

- there exists an internal predicted start site with score >= `cryptkeeper_fail_score`
- internal means start position:
  - >= `cryptkeeper_ignore_first_nt`
  - and < (sequence_length - `cryptkeeper_ignore_last_nt`)

The TSV includes:

- `cryptkeeper_internal_site_count`
- `cryptkeeper_max_internal_score`
- `cryptkeeper_internal_site_positions`

Counts are reported for debugging, but filtering is based on the maximum internal score, not the count.

## Performance: pooling (two-phase screening)

CryptKeeper can be slow, so codonopt uses a pooling strategy:

1. Generate a pool of candidates without CryptKeeper (fast).
2. Run CryptKeeper on the pool (slow).
3. Keep the first `n` candidates that pass.

Pool knobs:

- `cryptkeeper_pool_factor` (default 3.0)  
  Pool target is approximately `ceil(n * pool_factor)`.
- `cryptkeeper_max_pool` (default 200)  
  Hard cap to limit runtime.

If you request `n=10`, default pool target is about 30.  
If you request `n=100`, default pool target caps at 200 unless you raise `cryptkeeper_max_pool`.

## If CryptKeeper is too strict (not enough sequences passing)

Use one or more of the following adjustments:

1. Make it more permissive by raising the fail threshold:
   - Increase `cryptkeeper_fail_score` (e.g., 8000, 10000, 15000)

2. Consider fewer sites by increasing the reporting cutoff:
   - Increase `cryptkeeper_rbs_score_cutoff` (e.g., 1000)

3. Increase search coverage:
   - Increase `cryptkeeper_pool_factor` (e.g., 5.0, 10.0)
   - Increase `cryptkeeper_max_pool` (e.g., 500, 1000)
   - Increase `max_tries_per_replicate` if generation itself is failing

## If CryptKeeper is too permissive (you want stricter filtering)

1. Lower `cryptkeeper_fail_score` (e.g., 4000–5000)
2. Lower `cryptkeeper_rbs_score_cutoff` (e.g., 200–300) to consider more sites
3. Reduce ignore windows if you want to count near-end or near-start internal sites

## Quick calibration trick

If you don’t know what your score distribution looks like yet:

- Set `cryptkeeper_fail_score` extremely high (e.g., 1e9) so almost everything passes.
- Run once and inspect `cryptkeeper_max_internal_score` in `metrics.tsv`.
- Choose a fail threshold that removes the extreme tail (e.g., keep most sequences but reject the worst 5–10%).

---

## Batch table columns (common)

Required:

- `sequence`
- `codon_table`

Common optional columns:

- `codon_table_sheet`
- `n`
- `seed`
- `optimization_mode` (kleinbub|strict)
- `gc_min`, `gc_max`
- `avoid_codons` (use `|`)
- `avoid_motifs` (use `|`)
- `max_homopolymer`
- `min_codon_fraction`
- `max_tries_per_replicate`
- `kleinbub_search_limit`
- `backtrack_window`

CryptKeeper columns:

- `cryptkeeper_enable` (0/1)
- `cryptkeeper_rbs_score_cutoff`
- `cryptkeeper_fail_score`
- `cryptkeeper_ignore_first_nt`
- `cryptkeeper_ignore_last_nt`
- `cryptkeeper_threads`
- `cryptkeeper_pool_factor`
- `cryptkeeper_max_pool`

Example row:

    sequence,codon_table,codon_table_sheet,n,seed,optimization_mode,gc_min,gc_max,max_homopolymer,min_codon_fraction,cryptkeeper_enable,cryptkeeper_rbs_score_cutoff,cryptkeeper_fail_score,cryptkeeper_ignore_first_nt,cryptkeeper_ignore_last_nt,cryptkeeper_pool_factor,cryptkeeper_max_pool,cryptkeeper_threads
    /data/input.fasta,/data/Ecoli.xlsx,Ecoli1,10,42,kleinbub,0.40,0.60,5,0.05,1,500,6000,30,30,5.0,500,4

---

## Notes

- Maximum supported input FASTA records per run: 1000
- FASTA inputs in batch rows can contain multiple sequences
- Outputs are consolidated into a single FASTA + TSV
