# 🐸 FROGGER  
**Fragmented Recombined Orfs, Golden Gate Excision Ready**

FROGGER is a deterministic, Dockerized bioinformatics pipeline for:

1. Codon optimization and back-translation via `codonopt`
2. Splitting optimized ORFs into user-defined fragments
3. Adding Golden Gate cloning adapters to each fragment
4. Reassembling fragments into all valid recombinant ORFs while preserving fragment order
5. Validating the final post-ligation ORFs against biological and user-defined constraints

FROGGER is designed for library construction, modular gene design, and high-throughput cloning workflows where correctness, auditability, and reproducibility matter.

---

## Why FROGGER exists

Codon optimization alone is not sufficient for modern cloning pipelines.

Once ORFs are fragmented, recombined, and ligated, new biological artifacts can be introduced, including:

- forbidden motifs (for example restriction sites)
- unexpected homopolymers
- frameshifts or internal stop codons
- GC drift
- junction-created patterns that never existed in the original sequences

FROGGER explicitly models the post-ligation ORF, not just the individual fragments, and verifies that nothing unwanted was introduced during the design process.

---

## Core capabilities

### Codon optimization

- Uses codonopt version 1.2 (vendored as a git submodule)
- Supports:
  - protein FASTA or CDS FASTA
  - batch CSV or TSV input
  - deterministic runs via seed
  - bounded backtracking mode called `kleinbub` or strict optimization
  - CryptKeeper filtering as an optional screening step

### Fragmentation

- User-defined cut points, either global for all sequences or per-gene
- Optional frame enforcement so cut points must be multiples of 3
- Produces:
  - core coding fragments used for reassembly and post-ligation validation
  - cloning-ready fragments that include user-defined adapters

### Golden Gate preparation

- Adds user-specified adapters to each fragment
- Adapters are output-only and are not scanned or validated
- Clean separation between:
  - cloning fragments that include adapters
  - post-ligation ORFs that exclude adapters and represent the biological coding sequence

### Recombinant assembly

- Reassembles fragments into all valid combinations
- Fragment order is preserved
- Fragment sources may mix across genes
- Deterministic random sampling is used when the number of combinations is too large to enumerate exhaustively

### Final validation on the post-ligation ORF

- Scans the assembled post-ligation ORF for:
  - forbidden motifs
  - GC bounds
  - homopolymer limits
  - internal stop codons
  - frame integrity
- Reports exact motif locations and sequence context
- Ensures biological correctness of every accepted construct

---

## Repository structure

This repository uses two repos in one workspace via a git submodule.

```text
frogger-repo/
  codonopt/        # git submodule (codonopt v1.2)
  frogger/         # FROGGER pipeline
    Dockerfile
    requirements.txt
    FROGGER.schema.md
    frogger/
      cli.py
      pipeline.py
      util.py
      splitters.py
      gg_adapters.py
      recombine.py
      checks.py
      report.py
      io.py
```

---

## Installation and setup

### Clone (important: submodules)

```bash
git clone --recurse-submodules https://github.com/<your-org>/<frogger-repo>.git
cd <frogger-repo>
```

If you already cloned without submodules, run:

```bash
git submodule update --init --recursive
```

---

## Build the Docker image

From the repository root, run:

```bash
docker build -t frogger:0.1 -f frogger/Dockerfile .
```

This image bundles:

- FROGGER
- codonopt version 1.2
- all runtime dependencies

No conda. No Docker-in-Docker.

---

## Running FROGGER

### Basic invocation

```bash
docker run --rm \
  -v "$PWD/in:/in" \
  -v "$PWD/out:/out" \
  frogger:0.1 \
  run --config /in/pipeline.yaml --outdir /out
```

### Required inputs

- `pipeline.yaml` for full pipeline configuration
- input FASTA (protein or CDS)
- codon table XLSX (if not using batch mode)

---

## Configuration overview

The pipeline is fully controlled by a single YAML file.

Key sections:

```yaml
inputs:            # input FASTA
codonopt:          # codon optimization parameters
splits:            # fragment boundaries
golden_gate:       # cloning adapters
reassembly:        # combinatorics and sampling
final_checks:      # ORF validation
final_forbidden:   # motifs banned post-ligation
outputs:           # filenames
```

See `frogger/FROGGER.schema.md` for the complete schema.

---

## Outputs

All outputs are written to `--outdir`.

### 1) Codonopt outputs

Located in:

```text
out/codonopt_out/
  optimized_sequences.fasta
  metrics.tsv
```

### 2) Cloning fragments FASTA

Output file:

```text
fragments_with_gg.fasta
```

Each sequence:

- includes Golden Gate adapters
- is ready for synthesis or cloning
- header includes gene ID and fragment index

### 3) Reassembled ORFs FASTA

Output file:

```text
reassembled_orfs.fasta
```

Each sequence:

- represents the post-ligation coding sequence
- contains no adapters
- header includes construct ID, fragment provenance, and PASS or FAIL status

### 4) Main report

Output file:

```text
report.tsv
```

One row per construct, including:

- fragment provenance
- ORF length
- GC fraction
- homopolymer statistics
- pass or fail status
- failure reasons
- seeds used
- codonopt version

### 5) Motif hit report

Output file:

```text
junction_hits.tsv
```

Lists every forbidden motif detected, with:

- exact position
- sequence context
- source construct

---

## Determinism and reproducibility

FROGGER is fully deterministic when seeds are provided.

- Codon optimization seed
- Reassembly sampling seed
- Config hash (implicit via git)

Re-running with the same config, codonopt version, and seeds will produce bit-identical outputs.

---

## Common use cases

- Modular Golden Gate library construction
- Large-scale ORF recombination experiments
- Avoidance of restriction sites across ligation junctions
- Synthetic gene design under strict biological constraints
- Auditable, production-grade DNA design workflows

---

## Design philosophy

- Post-ligation biology matters
- Fragments are not genes
- Constraints must be re-checked after every transformation
- Sampling must be deterministic
- Pipelines must be auditable

FROGGER encodes these principles explicitly.

---

## Development notes

- `codonopt` is treated as read-only in this repo
- To update codonopt:

```bash
cd codonopt
git fetch
git checkout <new-tag-or-commit>
cd ..
git add codonopt
git commit -m "Update codonopt submodule"
```

---

## License and status

- FROGGER: project-defined license
- codonopt: see the codonopt repository
