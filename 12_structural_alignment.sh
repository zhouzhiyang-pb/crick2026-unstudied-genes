#!/bin/bash
# Step 12: structural alignment of every focused structure against the GOI,
# reporting TM-scores (Foldseek).
#
# Input : <GENE>/alphafold/*.pdb   (from step 11)
# Output: <GENE>/tm_scores.tsv     (query, target, TM-score per pair)
#
# Requires: foldseek (installed in step 05, into the bioinfo env).

set -e

if ! command -v foldseek >/dev/null 2>&1; then
    echo "foldseek not found -- install it first (step 05): conda install -n bioinfo foldseek -c bioconda"
    exit 1
fi

for gene in SPBC115.03 SPAPB1A10.12c; do
    fasta="${gene}/focused_sequences.fasta"
    alpha_dir="${gene}/alphafold"
    [ -f "$fasta" ] || { echo "missing input: $fasta"; continue; }
    [ -d "$alpha_dir" ] || { echo "missing directory: $alpha_dir (run step 11 first)"; continue; }

    # The GOI is the sequence whose label contains "our_gene".
    goi_ac=$(grep '^>.*our_gene' "$fasta" | head -1 | sed 's/^>.*|//')
    goi_pdb="${alpha_dir}/${goi_ac}.pdb"

    if [ ! -f "$goi_pdb" ]; then
        echo "missing GOI structure: $goi_pdb"
        continue
    fi

    out="${gene}/tm_scores.tsv"
    tmp="${gene}/foldseek_tmp"
    mkdir -p "$tmp"

    # Search the GOI against every focused structure and keep the TM-score.
    foldseek easy-search "$goi_pdb" "$alpha_dir" "$out" "$tmp" \
        --format-output query,target,alntmscore

    # Foldseek leaves intermediate database files in the temp dir; remove them.
    rm -rf "$tmp"

    echo "$gene: wrote $out"
done
