#!/usr/bin/env python3
"""
Step 13: generate the candidate-residue list for mapping onto the 3D structure
(the input to the manual ChimeraX step).

Combines three sources into one list of GOI residue positions, in UniProt
numbering (which matches the AlphaFold PDB):

  1. Conservation (step 10) -- the most-conserved, low-gap alignment columns,
     mapped to the GOI's residue number. Computed on the UNtrimmed MSA so the
     residue numbers match UniProt / the AlphaFold PDB.
  2. Known residues (step 03) -- UniProt-annotated residues, e.g. the covalent
     FAD-histidine or binding sites.
  3. Literature -- catalytic / binding residues reported in papers, mapped from
     the characterized homolog onto the GOI by pairwise sequence alignment.

Input (all inside <GENE>/):
  focused_sequences.fasta       (GOI + homolog sequences, unaligned)
  focused_sequences_msa.fasta   (untrimmed MSA, for conservation)
  annotations_table.csv         (step 03)
  residues_from_papers.tsv      (literature; columns: accession, residue, role, pmid)

Output:
  <GENE>/candidate_residues.csv  (goi_residue, aa, source, evidence)

Usage:
  python3 13_candidate_residues.py <GENE> [--top 20] [--max-gap 0.3]
"""

import argparse
import csv
import math
import re

from Bio import Align, AlignIO

GAP_CHARS = set("-.")
GOI_ACCESSION = {"SPBC115.03": "O42896", "SPAPB1A10.12c": "Q9HDX8"}

# BLOSUM62 background frequencies (used by the JSD conservation score).
BACKGROUND = {
    "A": 0.074, "R": 0.052, "N": 0.045, "D": 0.054, "C": 0.025,
    "Q": 0.034, "E": 0.054, "G": 0.074, "H": 0.026, "I": 0.068,
    "L": 0.099, "K": 0.058, "M": 0.025, "F": 0.047, "P": 0.039,
    "S": 0.057, "T": 0.051, "W": 0.013, "Y": 0.032, "V": 0.073,
}
_BG_ENTROPY = -sum(f * math.log2(f) for f in BACKGROUND.values())

# Global alignment, used to map a homolog's residue onto the GOI.
_aligner = Align.PairwiseAligner()


def entropy(probs):
    return -sum(p * math.log2(p) for p in probs if p > 0)


def column_jsd(column):
    """JSD conservation of one column; returns (score, most_common, gap_frac)."""
    counts = {}
    n_gap = 0
    for residue in column:
        residue = residue.upper()
        if residue in GAP_CHARS:
            n_gap += 1
            continue
        counts[residue] = counts.get(residue, 0) + 1
    total = sum(counts.values())
    if total == 0:
        return 0.0, "-", 1.0
    p = {aa: counts.get(aa, 0) / total for aa in BACKGROUND}
    mean = {aa: (p[aa] + BACKGROUND[aa]) / 2 for aa in BACKGROUND}
    jsd = entropy(mean.values()) - (entropy(p.values()) + _BG_ENTROPY) / 2
    gap_fraction = n_gap / len(column)
    return jsd * (1 - gap_fraction), max(counts, key=counts.get), gap_fraction


def read_fasta(path):
    """Return {accession: sequence} from a FASTA with '>label|accession' headers."""
    seqs = {}
    name = None
    parts = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(parts)
            header = line[1:]
            name = header.split("|")[-1] if "|" in header else header
            parts = []
        elif line:
            parts.append(line)
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def residue_number_at(goi_seq, col):
    """1-based residue number of the GOI at this alignment column (None if gap)."""
    if goi_seq[col] in GAP_CHARS:
        return None
    return sum(1 for i in range(col + 1) if goi_seq[i] not in GAP_CHARS)


def residue_at(seq, pos):
    """Return the amino acid at UniProt position pos (1-based)."""
    n = 0
    for c in seq:
        if c not in GAP_CHARS:
            n += 1
            if n == pos:
                return c
    return "?"


def map_position(homolog_seq, goi_seq, pos):
    """
    Map a homolog residue position (1-based) to the GOI's residue number
    (1-based) using a global pairwise alignment. Returns None if the homolog
    residue aligns to a gap in the GOI.
    """
    if pos < 1 or pos > len(homolog_seq):
        return None
    aln = _aligner.align(homolog_seq, goi_seq)[0]
    h, g = aln[0], aln[1]  # aligned homolog / GOI (with gaps)
    h_count = 0
    for col in range(len(h)):
        if h[col] in GAP_CHARS:
            continue
        h_count += 1
        if h_count == pos:
            if g[col] in GAP_CHARS:
                return None
            return sum(1 for i in range(col + 1) if g[i] not in GAP_CHARS)
    return None


def read_literature(tsv_path):
    """Read residues_from_papers.tsv -> list of (accession, [positions], role, pmid)."""
    entries = []
    for line in open(tsv_path):
        line = line.rstrip("\n")
        if not line or line.startswith("#") or line.startswith("accession"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        accession, residue, role, pmid = parts[0], parts[1], parts[2], parts[3]
        positions = []
        for token in residue.split(","):
            token = token.strip()
            if "-" in token:  # a range like "75-96"
                a, b = token.split("-")
                positions.extend(range(int(a), int(b) + 1))
            elif token.isdigit():
                positions.append(int(token))
        entries.append((accession, positions, role, pmid))
    return entries


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gene", choices=["SPBC115.03", "SPAPB1A10.12c"])
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--max-gap", type=float, default=0.3)
    args = ap.parse_args()
    gene = args.gene
    goi_ac = GOI_ACCESSION[gene]

    # Ungapped sequences (step 04 fasta), for literature mapping + residue identity.
    seqs = read_fasta(f"{gene}/focused_sequences.fasta")
    goi_ungapped = seqs[goi_ac]

    candidates = []  # (goi_residue, aa, source, evidence)

    # --- 1. conservation (on the untrimmed MSA) ---
    aln = AlignIO.read(f"{gene}/focused_sequences_msa.fasta", "fasta")
    msa_goi = next(str(rec.seq) for rec in aln if "our_gene" in rec.id)
    cons = []
    for col in range(aln.get_alignment_length()):
        score, _, gap = column_jsd(aln[:, col])
        if gap > args.max_gap:
            continue
        res = residue_number_at(msa_goi, col)
        if res is not None:
            cons.append((score, res, msa_goi[col]))
    cons.sort(key=lambda x: -x[0])
    for score, res, aa in cons[: args.top]:
        candidates.append([res, aa, "conservation", f"JSD={score:.3f}"])

    # --- 2. known residues (annotations, already in UniProt numbering) ---
    for r in csv.reader(open(f"{gene}/annotations_table.csv")):
        if r and r[0] == goi_ac:
            for feature in (r[7] or "").split(";"):
                feature = feature.strip()
                m = re.search(r"(?:residue|site)\s+(\d+)", feature)
                if m:
                    pos = int(m.group(1))
                    candidates.append([pos, residue_at(goi_ungapped, pos), "annotation", feature])

    # --- 3. literature residues (mapped from homologs) ---
    for accession, positions, role, pmid in read_literature(f"{gene}/residues_from_papers.tsv"):
        homolog_seq = seqs.get(accession)
        if homolog_seq is None:
            print(f"  warning: {accession} not in focused_sequences.fasta; skipped")
            continue
        for pos in positions:
            mapped = map_position(homolog_seq, goi_ungapped, pos)
            if mapped is None:
                continue
            candidates.append([mapped, residue_at(goi_ungapped, mapped), "literature",
                               f"{role} (from {accession} #{pos}, PMID {pmid})"])

    # --- merge, dedup by residue number; keep the most-informative source ---
    # annotation (known UniProt residue) > literature (from papers) > conservation (predicted)
    PRIORITY = {"annotation": 0, "literature": 1, "conservation": 2}
    best = {}
    for res, aa, src, ev in candidates:
        if res not in best or PRIORITY[src] < PRIORITY[best[res][2]]:
            best[res] = [res, aa, src, ev]
    out_rows = list(best.values())

    out = f"{gene}/candidate_residues.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["goi_residue", "aa", "source", "evidence"])
        w.writerows(sorted(out_rows, key=lambda x: int(x[0])))

    n_lit = sum(1 for r in out_rows if r[2] == "literature")
    print(f"wrote {len(out_rows)} candidate residues -> {out}")
    print(f"  (conservation + annotation + {n_lit} literature-mapped)")


if __name__ == "__main__":
    main()
