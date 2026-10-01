import re
import math
from collections import Counter
from itertools import combinations

import numpy as np
import pandas as pd
from Bio import SeqIO

import config

VALID_BASES = ("A", "C", "G", "T")
HAPLOGROUP_QUALITY_THRESHOLD = 0.85


def load_alignment(fasta_path=config.OUTPUT_FASTA):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    seqs = [str(r.seq).upper() for r in records]
    ids = [r.id for r in records]

    lengths = Counter(len(s) for s in seqs)
    consensus_len = lengths.most_common(1)[0][0]
    keep = [(i, s) for i, s in zip(ids, seqs) if len(s) == consensus_len]
    dropped = len(seqs) - len(keep)
    if dropped:
        print(f"Dropped {dropped} sequence(s) with inconsistent alignment length")

    ids = [i for i, _ in keep]
    seqs = [s for _, s in keep]
    matrix = np.array([list(s) for s in seqs])
    return ids, seqs, matrix


def column_base_counts(matrix):
    counts = {b: (matrix == b).sum(axis=0) for b in VALID_BASES}
    n_valid = sum(counts.values())
    return counts, n_valid


def segregating_sites(matrix):
    counts, n_valid = column_base_counts(matrix)
    total_pairs = n_valid * (n_valid - 1) / 2
    same_pairs = sum(c * (c - 1) / 2 for c in counts.values())
    diff_pairs = total_pairs - same_pairs
    seg_mask = diff_pairs > 0
    S = int(seg_mask.sum())
    return S, diff_pairs, total_pairs, seg_mask


class DisjointSet:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx != ry:
            self.parent[rx] = ry


def sequences_compatible(a, b):
    informative = (a != "N") & (b != "N")
    if not informative.any():
        return True
    return bool(np.all(a[informative] == b[informative]))


def haplotype_stats(matrix, seg_mask):
    n = matrix.shape[0]
    sub = matrix[:, seg_mask]

    dsu = DisjointSet(n)
    for i, k in combinations(range(n), 2):
        if sequences_compatible(sub[i], sub[k]):
            dsu.union(i, k)

    groups = Counter(dsu.find(i) for i in range(n))
    H = len(groups)
    probs = np.array([c / n for c in groups.values()])
    j = float(np.sum(probs ** 2))
    Hd = (n / (n - 1)) * (1 - j) if n > 1 else 0.0
    return H, Hd, j


def nucleotide_diversity(diff_pairs, total_pairs, seq_length):
    valid = total_pairs > 0
    site_pi = np.zeros_like(diff_pairs, dtype=float)
    site_pi[valid] = diff_pairs[valid] / total_pairs[valid]

    L_used = int(valid.sum())
    k_hat_raw = float(site_pi[valid].sum())
    pi_per_site = k_hat_raw / L_used if L_used > 0 else 0.0
    return pi_per_site, k_hat_raw, L_used


def wattersons_theta(S, n, seq_length):
    a1 = sum(1 / i for i in range(1, n))
    theta_w_raw = S / a1 if a1 > 0 else 0.0
    theta_w_per_site = theta_w_raw / seq_length
    return theta_w_per_site, theta_w_raw, a1


def tajimas_d(S, k_hat_raw, theta_w_raw, n, a1):
    if S == 0 or n < 3:
        return 0.0

    a2 = sum(1 / (i ** 2) for i in range(1, n))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n ** 2 + n + 2) / (9 * n * (n - 1))
    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 ** 2)
    e1 = c1 / a1
    e2 = c2 / (a1 ** 2 + a2)

    variance = e1 * S + e2 * S * (S - 1)
    if variance <= 0:
        return 0.0

    D = (k_hat_raw - theta_w_raw) / math.sqrt(variance)
    return D


def major_haplogroup(label):
    if not label or label.upper() in ("NA", "N/A", ""):
        return "Unclassified"
    match = re.match(r"^([A-Za-z]+)(\d*)", label.strip())
    if not match:
        return label
    letters, digits = match.group(1), match.group(2)
    if letters.upper() == "L" and digits:
        return f"{letters}{digits[0]}"
    return letters


def haplogroup_breakdown(haplogroups_path=None, quality_threshold=HAPLOGROUP_QUALITY_THRESHOLD):
    path = haplogroups_path or getattr(config, "HAPLOGREP_OUTPUT", config.WORKDIR / "haplogroups.txt")
    if not path.exists():
        print(f"No haplogroups file found at {path}, skipping haplogroup breakdown")
        return pd.Series(dtype=int)

    df = pd.read_csv(path, sep="\t")
    hg_col = next((c for c in ("Haplogroup",) if c in df.columns), df.columns[-1])
    quality_col = next((c for c in ("Quality",) if c in df.columns), None)

    if quality_col is None:
        print("No 'Quality' column found in haplogroups.txt — skipping confidence filtering, "
              "check the actual column header if this is unexpected")

    def classify(row):
        if quality_col is not None:
            try:
                quality = float(row[quality_col])
                if quality < quality_threshold:
                    return "Low_Confidence_Call"
            except (TypeError, ValueError):
                pass
        return major_haplogroup(row[hg_col])

    majors = df.apply(classify, axis=1)
    return majors.value_counts()


def run():
    ids, seqs, matrix = load_alignment()
    n = len(seqs)
    seq_length = matrix.shape[1]

    S, diff_pairs, total_pairs, seg_mask = segregating_sites(matrix)
    H, Hd, j = haplotype_stats(matrix, seg_mask)
    pi_per_site, k_hat_raw, L_used = nucleotide_diversity(diff_pairs, total_pairs, seq_length)
    theta_w_per_site, theta_w_raw, a1 = wattersons_theta(S, n, seq_length)
    D = tajimas_d(S, k_hat_raw, theta_w_raw, n, a1)
    hg_counts = haplogroup_breakdown()

    rows = [
        {"Metric": "N (sequences)", "Value": n},
        {"Metric": "Alignment length (L)", "Value": seq_length},
        {"Metric": "Sites used for Pi (L_used)", "Value": L_used},
        {"Metric": "S (segregating sites)", "Value": S},
        {"Metric": "H (unique haplotypes, N-aware)", "Value": H},
        {"Metric": "Hd (haplotype diversity)", "Value": round(Hd, 6)},
        {"Metric": "j (homozygosity)", "Value": round(j, 6)},
        {"Metric": "Pi (nucleotide diversity, per site)", "Value": round(pi_per_site, 6)},
        {"Metric": "Theta_W (Watterson, per site)", "Value": round(theta_w_per_site, 6)},
        {"Metric": "Tajima's D", "Value": round(D, 6)},
    ]

    for hg, count in hg_counts.sort_values(ascending=False).items():
        rows.append({"Metric": f"Haplogroup_{hg}", "Value": int(count)})

    out_path = config.WORKDIR / "popgen_summary.csv"
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    run()