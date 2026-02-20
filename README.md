# 🐸 FROGGER  
**Fragmented Recombined Orfs, Golden Gate Excision Ready**

FROGGER is a unified, Dockerized **Golden Gate fragment library generator** for protein engineering workflows. It provides a single CLI with two complementary modes:

- **`frogger shuffle`** — recombination / shuffling across homologous proteins (legacy FROGGER)
- **`frogger mutate`** — site-saturation / targeted mutagenesis with **single-mutation** fragment pooling + assembly plan (ported from SLOGGER)

Both workflows converge on the same shared Golden Gate (GG) packaging logic: junction-aware overhang selection, forced codons to realize junction overhangs, adapter/primer construction, and audit reporting.

---

## What FROGGER does

Given protein sequences (and a configuration describing split points and GG constraints), FROGGER produces GG-ready fragment libraries that:

- are **codon-optimized** (via `codonopt`)
- are split into **Golden Gate fragments** at user-defined AA cut points
- use **junction-aware overhang selection** with a crosstalk scoring matrix + filters + beam search
- **force junction codons** so chosen overhangs are realizable in the final DNA
- generate **GG cloning-ready fragments** with enzyme sites, spacers, overhangs, and primers
- (mutate) produce **fragment pools** that assemble into **only single-mutant constructs**
- emit audit artifacts (overhang choice, crosstalk, reports)

---

## Quick start (Docker)

### Build

From the repo root:

    docker build -t frogger:1.1.1 .

### Run (recommended workflow)

Mount a local folder that contains your config and input FASTA files:

    docker run --rm -v "/path/to/examples:/examples" frogger:1.0 shuffle --config /examples/config_unified.yaml --outdir /examples/out_shuffle

    docker run --rm -v "/path/to/examples:/examples" frogger:1.0 mutate --config /examples/config_unified.yaml --outdir /examples/out_mutate

Notes:
- Use absolute paths for the `-v` mount.
- Output directories should be distinct per run to avoid overwrites.

---

## Installation (from source)

FROGGER is primarily used via Docker, but can be run locally if you prefer.

### Clone with submodules (codonopt)

This repo expects `codonopt` at the repo root path `./codonopt` (as a git submodule).

Clone with submodules:

    git clone --recurse-submodules <YOUR_FROGGER_REPO_URL>
    cd frogger

Or initialize after cloning:

    git submodule update --init --recursive

### Python environment

From the repo root:

    python -m venv .venv
    source .venv/bin/activate
    pip install -r frogger_app/requirements.txt

Then run the CLI (see “Usage”).

---

## Repository layout

- `frogger_app/`
  - `frogger/` — main Python package
    - `modules/shuffle/` — shuffle workflow
    - `modules/mutate/` — mutate workflow
    - `common/gg/` — shared Golden Gate helpers (splitters, packager)
    - `primers.py`, `overhangs.py`, `junction_overhangs.py`, `gg_adapters.py` — core shared GG logic
- `codonopt/` — codon optimization engine (git submodule)
- `examples/` — example configs and inputs

---

## Unified configuration

A single YAML config (commonly `config_unified.yaml`) drives either mode. Typical blocks:

- `inputs` — FASTA paths and input interpretation
- `codonopt` — codon optimization settings
- `splits` — AA cut points (shuffle uses directly; mutate uses to derive nt cuts)
- `golden_gate` — overhang matrix path, filters, beam search settings
  - `golden_gate.fragments` — enzyme, spacers, etc.
  - `golden_gate.max_worst_crosstalk` — worst allowed crosstalk threshold (default often `0.1`)
- `primers` — fixed or generated; motif avoidance / kmer avoidance
- `outputs` — output filenames (TSVs/FASTAs)
- `repair` — optional post-generation QC + repair (see below)

---

## Usage

FROGGER exposes a single `frogger` command with subcommands.

### Shuffle (recombination / shuffling)

    frogger shuffle --config /path/to/config_unified.yaml --outdir /path/to/out_shuffle

Shuffle produces GG-ready fragments representing recombined ORFs across a set of homologous proteins, with overhang selection optimized to minimize crosstalk.

### Mutate (site-saturation / targeted mutagenesis)

    frogger mutate --config /path/to/config_unified.yaml --outdir /path/to/out_mutate

Mutate produces:
- wild-type fragments
- mutant fragments
- fragment pools such that assembly yields only **single-mutation** constructs
- an `assembly_plan.tsv` mapping each variant to its required fragment combination

Important behavior:
- Full saturation (default) skips mutating amino acid position 1 (start codon) by default.

---

## Golden Gate packaging (shared behavior)

Both workflows share the same GG packaging steps:

1. Enumerate overhang candidates from synonymous codon choices at junction AAs
2. Optimize a global overhang assignment using:
   - pairing/crosstalk matrix scoring
   - filters
   - beam search
3. Force codons at junction AA indices to realize chosen overhangs
4. Build cloning sequences:
   - enzyme recognition sites
   - overhangs
   - spacers
   - primers (fixed or generated)
5. Emit audit files, including:
   - selected overhangs per junction
   - pairwise crosstalk scores
   - fail-fast if worst crosstalk exceeds `golden_gate.max_worst_crosstalk`

---

## Post-generation QC + repair (optional)

Some outputs (especially shuffle) can occasionally violate downstream synthesis/assembly constraints such as:
- excessive homopolymer runs
- forbidden motifs
- internal restriction enzyme sites

FROGGER supports an optional **QC + repair** pass that:
1. Scans fragment cores (and/or ORFs, depending on implementation) for constraint violations
2. If failures exist, re-runs codon optimization **only on failing sequences**
3. Validates **AA preservation** (translation must match exactly)
4. Re-applies junction forced codons (when repairing fragment cores)
5. Regenerates downstream packaging so results remain Golden Gate compatible
6. Writes a clear `repair_report.tsv` describing what changed and what still fails

### Example config block

    repair:
      enable: true
      max_rounds: 1
      qc:
        max_homopolymer: 5
        forbidden_motifs: ["AAAAAA", "TTTTTT"]
        internal_sites:
          - {name: "BsaI", site: "GGTCTC"}
          - {name: "BsmBI", site: "CGTCTC"}

### Repair report

When enabled and any repairs are attempted, FROGGER writes:

- `repair_report.tsv` (name configurable via `outputs.repair_report_tsv`)

Typical columns include:
- fragment/record id
- round number
- failures before/after
- whether AA was preserved
- whether the nucleotide sequence changed
- codonopt version (if available)

---

## Outputs (typical)

Exact filenames are controlled by the `outputs` section in the config, but commonly include:

- FASTA files:
  - GG-ready fragments (with adapters/primers)
  - reconstructed ORFs (where applicable)
  - mutant variant sequences (mutate)
- TSV reports:
  - pipeline report summary
  - junction/overhang audit
  - overhang crosstalk table
  - (mutate) assembly plan and pooling tables
  - (optional) repair report

---

## Development notes

- `codonopt` is a git submodule at `./codonopt`.
- If you update the codonopt submodule pointer:

    cd codonopt
    git fetch
    git checkout <new_commit_or_tag>
    cd ..
    git add codonopt
    git commit -m "Bump codonopt submodule"
    git push

- For reproducibility, keep deterministic seeds consistent in config and/or pipeline defaults.

---

## Troubleshooting

### Submodule not present
If builds fail because `codonopt` is missing, initialize submodules:

    git submodule update --init --recursive

### Docker build is missing codonopt
If you add a `.dockerignore`, ensure it does not exclude `codonopt/`.

---
