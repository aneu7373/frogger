from pathlib import Path
import pandas as pd
from Bio import SeqIO


def ensure_outdir(path: str) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def read_fasta_records(path: Path):
    return [(r.id, str(r.seq)) for r in SeqIO.parse(str(path), "fasta")]


def write_fasta_records(path: Path, records):
    with open(path, "w", encoding="utf-8") as f:
        for h, s in records:
            f.write(f">{h}\n")
            for i in range(0, len(s), 80):
                f.write(s[i : i + 80] + "\n")


def write_tsv(path: Path, rows):
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)
