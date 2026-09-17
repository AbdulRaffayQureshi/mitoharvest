import os
import csv
import hashlib
import config

FIELDNAMES = ["accession", "region", "length", "pubmed_id", "pubmed_title", "raw_sequence", "aligned_sequence", "description"]


def sequence_hash(seq: str) -> str:
    return hashlib.md5(seq.upper().encode()).hexdigest()


def load_checkpoint(csv_path=config.OUTPUT_CSV):
    seen_accessions = set()
    seen_hashes = set()
    count = 0

    if not os.path.exists(csv_path):
        return seen_accessions, seen_hashes, count

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            seen_accessions.add(row["accession"])
            seen_hashes.add(sequence_hash(row["raw_sequence"]))
            count += 1

    return seen_accessions, seen_hashes, count


class ResultWriter:
    def __init__(self, csv_path=config.OUTPUT_CSV, fasta_path=config.OUTPUT_FASTA, resume=False):
        csv_exists = resume and os.path.exists(csv_path)
        fasta_exists = resume and os.path.exists(fasta_path)

        self.csv_file = open(csv_path, "a" if csv_exists else "w", newline="")
        self.csv_writer = csv.DictWriter(self.csv_file, fieldnames=FIELDNAMES)
        if not csv_exists:
            self.csv_writer.writeheader()
        self._flush(self.csv_file)

        self.fasta_file = open(fasta_path, "a" if fasta_exists else "w")

    def _flush(self, f):
        f.flush()
        os.fsync(f.fileno())

    def write(self, record):
        self.csv_writer.writerow({k: record[k] for k in FIELDNAMES})
        self._flush(self.csv_file)
        self.fasta_file.write(f">{record['accession']} {record['region']}\n{record['aligned_sequence']}\n")
        self._flush(self.fasta_file)

    def close(self):
        self.csv_file.close()
        self.fasta_file.close()


