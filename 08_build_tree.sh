#!/bin/bash
# Step 8: build a maximum-likelihood tree from one or more trimmed alignments (IQ-TREE 3).
#
# Usage: bash 08_build_tree.sh <trimmed_msa.fasta> [<trimmed2.fasta> ...]
#   e.g. bash 08_build_tree.sh SPBC115.03/focused_sequences_msa_trimmed.fasta
#
# For each input <dir>/<X>.fasta, writes IQ-TREE files to <dir>/tree/
# (the main tree is <dir>/tree/<X>.treefile).
#
# Requires: iqtree3 (installed in step 05, Homebrew).

set -e

if ! command -v iqtree3 >/dev/null 2>&1; then
    echo "iqtree3 not found -- install it (step 05) and activate the env: conda activate bioinfo"
    exit 1
fi

if [ "$#" -eq 0 ]; then
    echo "usage: bash 08_build_tree.sh <trimmed_msa.fasta> [more inputs...]"
    exit 1
fi

build_tree () {
    local input=$1
    local dir
    dir=$(dirname "$input")
    local base
    base=$(basename "${input%.fasta}")
    mkdir -p "$dir/tree"

    # -m MFP   : best amino-acid substitution model (ModelFinder Plus)
    # -bb 1000 : ultrafast bootstrap, 1000 replicates
    # -nt AUTO : use all CPU cores (needs the conda multicore build, not Homebrew's single-core)
    # -redo    : overwrite a previous run in the same folder
    iqtree3 -s "$input" -m MFP -bb 1000 -nt AUTO -pre "$dir/tree/$base" -redo

    echo "tree written to $dir/tree/$base.treefile"
}

for input in "$@"; do
    if [ ! -f "$input" ]; then
        echo "missing input file: $input"
        continue
    fi
    build_tree "$input"
done
