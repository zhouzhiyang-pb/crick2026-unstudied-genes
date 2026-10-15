#!/usr/bin/env python3
"""
Step 6: pairwise sequence identity (% ID) of every focused sequence against
the gene of interest (GOI).

The GOI is the sequence whose FASTA header label contains "our_gene".

Input : <GENE>/focused_sequences.fasta   (unaligned, from step 04)
Output: <GENE>/pairwise_percent_id.csv    (one row per focused sequence,
                                           sorted by % ID, highest first)

% ID is computed from a global pairwise alignment (Needleman-Wunsch, via
Biopython's PairwiseAligner):
    % ID = identical positions / alignment length * 100

Requires: biopython (already in the "bioinfo" conda env, as a clipkit dependency).

Run with the bioinfo env active:
  conda activate bioinfo
  python3 06_pairwise_percent_id.py
"""

import csv
import os

from Bio import Align

GENES = ["SPBC115.03", "SPAPB1A10.12c"]

# Global alignment (the default): align the two full-length sequences end to end.
_aligner = Align.PairwiseAligner()


def read_fasta(path):
    """Return a list of (label, accession, sequence) from a FASTA file."""
    records = []
    label = accession = None
    seq = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if label is not None:
                records.append((label, accession, "".join(seq)))
            header = line[1:]
            if "|" in header:
                label, accession = header.split("|", 1)
            else:
                label, accession = header, ""
            seq = []
        elif line:
            seq.append(line)
    if label is not None:
        records.append((label, accession, "".join(seq)))
    return records


def percent_identity(seq1, seq2):
    """Global % identity between two sequences (identical / aligned length)."""
    aln = _aligner.align(seq1, seq2)[0]  # best alignment
    a, b = aln[0], aln[1]
    matches = sum(1 for x, y in zip(a, b) if x == y)
    return matches / len(a) * 100.0


def main():
    for gene in GENES:
        fasta = f"{gene}/focused_sequences.fasta"
        if not os.path.isfile(fasta):
            print(f"missing input file: {fasta}")
            continue

        records = read_fasta(fasta)

        # The GOI is our own gene; its label contains "our_gene".
        goi = next((r for r in records if "our_gene" in r[0]), None)
        if goi is None:
            print(f"{gene}: no sequence labelled 'our_gene' found")
            continue

        goi_label, goi_ac, goi_seq = goi
        rows = []
        for label, accession, seq in records:
            if label == goi_label:
                continue
            rows.append([label, accession, f"{percent_identity(goi_seq, seq):.1f}"])

        rows.sort(key=lambda r: float(r[2]), reverse=True)

        out = f"{gene}/pairwise_percent_id.csv"
        with open(out, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["label", "accession", "percent_identity_vs_goi"])
            writer.writerows(rows)
        print(f"{gene}: wrote {len(rows)} rows -> {out}  (GOI = {goi_label} / {goi_ac})")


if __name__ == "__main__":
    main()
