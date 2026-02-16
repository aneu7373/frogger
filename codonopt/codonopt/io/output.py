import os

def write_sequence(seq, seq_id, out_dir, job_index):
    job_dir = os.path.join(out_dir, f"job_{job_index:04d}")
    os.makedirs(job_dir, exist_ok=True)
    out_file = os.path.join(job_dir, "optimized.fasta")
    with open(out_file, "w") as f:
        f.write(f">{seq_id}\n{seq}\n")
