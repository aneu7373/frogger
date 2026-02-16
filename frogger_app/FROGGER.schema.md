FROGGER pipeline config (YAML)

Top-level keys:

pipeline_version: "0.1"

inputs:
  fasta: "/in/inputs.fasta"
  input_type: "protein" | "cds"   # informational; codonopt will accept either

codonopt:
  # Inputs: choose ONE
  sequence: "/in/inputs.fasta"    # passed to --sequence
  batch_table: null               # passed to --batch-table (if provided, overrides many per-row)

  # Single-sequence mode codon table
  codon_table: "/in/codon_table.xlsx"
  codon_table_sheet: null

  out_subdir: "codonopt_out"
  verbose: false
  log_file: null

  avoid_codons: []                # list -> joined by '|'
  avoid_motifs: []                # list -> joined by '|'
  gc_min: 0.35
  gc_max: 0.65
  max_homopolymer: 5
  min_codon_fraction: 0.05
  n: 1
  seed: null

  optimization:
    mode: "kleinbub" | "strict"
    max_tries_per_replicate: 25
    kleinbub_search_limit: 200000
    backtrack_window: 10

  cryptkeeper:
    enable: false
    rbs_score_cutoff: 500.0
    threads: 1
    fail_score: 6000.0
    ignore_first_nt: 30
    ignore_last_nt: 30
    pool_factor: 3.0
    max_pool: 200

splits:
  global_cut_points: [300, 600]     # 0-based cut points into CDS
  per_gene_cut_points: {}          # optional gene_id -> [cut points]
  enforce_frame: true              # require cut points % 3 == 0

golden_gate:
  fragments:                        # one entry per fragment index (k entries)
    - left_adapter: "..."
      right_adapter: "..."
    - left_adapter: "..."
      right_adapter: "..."
    - left_adapter: "..."
      right_adapter: "..."

reassembly:
  max_constructs: 100000
  sample_n: 20000
  seed: null
  deduplicate: true

final_checks:
  junction_window: 24
  require_no_internal_stops: true
  require_length_multiple_of_3: true

final_forbidden:
  motifs: []                        # optional additional motifs to scan post-ligation

defaults:
  default_seed: 1337

outputs:
  fragments_fasta: "fragments_with_gg.fasta"
  reassembled_fasta: "reassembled_orfs.fasta"
  report_tsv: "report.tsv"
  junction_hits_tsv: "junction_hits.tsv"
