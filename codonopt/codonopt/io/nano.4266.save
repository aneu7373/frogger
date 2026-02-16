import pandas as pd

def write_outputs(results, fasta_path, tsv_path):
    records = []

    with open(fasta_path, "w") as fasta:
        for i, (seq, metrics) in enumerate(results, 1):
            fasta.write(f">variant_{i} score={metrics['score']:.3f}\n")
            fasta.write(seq + "\n")

            record = {"variant": i, **metrics}
            records.append(record)

    pd.DataFrame(records).to_csv(tsv_path, sep="\t", index=False)
