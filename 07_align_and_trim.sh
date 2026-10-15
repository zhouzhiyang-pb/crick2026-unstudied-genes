#!/bin/bash
# Step 7: align one or more sets of sequences (mafft) and trim each (ClipKIT).
#
# Usage: bash 07_align_and_trim.sh <unaligned.fasta> [<unaligned2.fasta> ...]
#   e.g. bash 07_align_and_trim.sh SPBC115.03/focused_sequences.fasta \
#                                SPAPB1A10.12c/focused_sequences.fasta
#
# For each input <X>.fasta, writes next to it:
#   <X>_msa.fasta          (aligned)
#   <X>_msa_trimmed.fasta  (trimmed)
#
# Requires:
#   mafft    -- system-wide (step 05, Homebrew)
#   clipkit  -- conda env "bioinfo" (step 05)
# Run with the bioinfo env active:
#   conda activate bioinfo
#   bash 07_align_and_trim.sh <input.fasta> ...

set -e

if ! command -v clipkit >/dev/null 2>&1; then
    echo "clipkit not found -- activate the environment first: conda activate bioinfo"
    exit 1
fi

if [ "$#" -eq 0 ]; then
    echo "usage: bash 07_align_and_trim.sh <unaligned.fasta> [more inputs...]"
    exit 1
fi

align_and_trim () {
    local input=$1
    local base="${input%.fasta}"
    local msa="${base}_msa.fasta"
    local trimmed="${base}_msa_trimmed.fasta"

    # Align with mafft.
    # --localpair uses L-INS-i, the most accurate strategy for sequence sets of this size.
    # --maxiterate 1000 maximizes the refinement rounds.
    mafft --localpair --maxiterate 1000 --quiet "$input" > "$msa"

    # Trim with ClipKIT. smart-gap removes gappy columns but keeps informative sites.
    clipkit "$msa" -m smart-gap -o "$trimmed"

    echo "aligned: $msa"
    echo "trimmed: $trimmed"
}

for input in "$@"; do
    if [ ! -f "$input" ]; then
        echo "missing input file: $input"
        continue
    fi
    align_and_trim "$input"
done
