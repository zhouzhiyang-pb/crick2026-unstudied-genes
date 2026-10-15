#!/usr/bin/env python3
"""
Step 10: per-column conservation score of an alignment.

Two methods:
  --method jsd      (default) Jensen-Shannon divergence, computed locally.
                    Tree-free: compares each column's amino-acid distribution
                    against the BLOSUM62 background (Capra & Singh 2007).
  --method consurf  prints instructions to upload the alignment to the ConSurf
                    web server (https://consurf.tau.ac.il), which computes the
                    Rate4Site (Bayesian) score and maps it onto a structure.

Usage:
  python3 10_conservation.py <msa.fasta> --method jsd
  python3 10_conservation.py <msa.fasta> --method consurf

Output (jsd): <msa>_conservation.csv
  columns: column, jsd_conservation, most_common_residue, gap_fraction

Run with the bioinfo env active:
  conda activate bioinfo
  python3 10_conservation.py <msa.fasta>
"""

import argparse
import csv
import math
import os
import sys

from Bio import AlignIO

GAP_CHARS = set("-.")

# Background amino-acid frequencies (BLOSUM62). JSD compares each column's
# distribution against this, so a rare residue (e.g. W) only counts as
# "conserved" if it occurs more often than expected by chance.
BACKGROUND = {
    "A": 0.074, "R": 0.052, "N": 0.045, "D": 0.054, "C": 0.025,
    "Q": 0.034, "E": 0.054, "G": 0.074, "H": 0.026, "I": 0.068,
    "L": 0.099, "K": 0.058, "M": 0.025, "F": 0.047, "P": 0.039,
    "S": 0.057, "T": 0.051, "W": 0.013, "Y": 0.032, "V": 0.073,
}
_BG_ENTROPY = -sum(f * math.log2(f) for f in BACKGROUND.values())


def entropy(probs):
    """Shannon entropy of a list of probabilities."""
    return -sum(p * math.log2(p) for p in probs if p > 0)


def column_conservation(column):
    """Return (JSD conservation, most common residue, gap fraction) for a column."""
    counts = {}
    n_gap = 0
    for residue in column:
        residue = residue.upper()
        if residue in GAP_CHARS:
            n_gap += 1
            continue
        counts[residue] = counts.get(residue, 0) + 1
    total = sum(counts.values())

    if total == 0:  # all-gap column
        return 0.0, "-", 1.0

    p = {aa: counts.get(aa, 0) / total for aa in BACKGROUND}
    mean = {aa: (p[aa] + BACKGROUND[aa]) / 2 for aa in BACKGROUND}
    jsd = entropy(mean.values()) - (entropy(p.values()) + _BG_ENTROPY) / 2

    gap_fraction = n_gap / len(column)
    score = jsd * (1 - gap_fraction)  # downweight gappy columns

    return score, max(counts, key=counts.get), gap_fraction


def run_consurf(msa_path):
    """Print instructions to run the ConSurf web tool (can't be done locally)."""
    print("ConSurf is a web tool; run it in the browser:")
    print("  1. go to https://consurf.tau.ac.il")
    print(f"  2. upload: {msa_path}")
    print("  3. choose the 'Bayesian' (Rate4Site) method")
    print("  4. download the per-residue conservation scores (1-9 scale)")
    print("  5. optionally map them onto a PDB / AlphaFold structure")


def run_jsd(msa_path):
    """Compute JSD conservation and write a CSV next to the input."""
    alignment = AlignIO.read(msa_path, "fasta")
    rows = []
    for i in range(alignment.get_alignment_length()):
        score, most_common, gap_frac = column_conservation(alignment[:, i])
        rows.append([i + 1, f"{score:.4f}", most_common, f"{gap_frac:.3f}"])

    out = os.path.splitext(msa_path)[0] + "_conservation.csv"
    with open(out, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["column", "jsd_conservation", "most_common_residue", "gap_fraction"])
        writer.writerows(rows)
    print(f"wrote {len(rows)} columns -> {out}")


def main():
    parser = argparse.ArgumentParser(description="Per-column conservation score (JSD or ConSurf).")
    parser.add_argument("input", help="path to the MSA in FASTA format")
    parser.add_argument("--method", choices=["jsd", "consurf"], default="jsd",
                        help="jsd = local Jensen-Shannon; consurf = ConSurf web instructions")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"missing input file: {args.input}")
        sys.exit(1)

    if args.method == "consurf":
        run_consurf(args.input)
    else:
        run_jsd(args.input)


if __name__ == "__main__":
    main()
