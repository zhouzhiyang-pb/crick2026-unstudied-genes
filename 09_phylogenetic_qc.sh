#!/bin/bash
# Step 9: phylogenetic signal QC (phykit).
#
# Usage: bash 09_phylogenetic_qc.sh <trimmed_msa.fasta> <tree.treefile>
#   e.g. bash 09_phylogenetic_qc.sh SPBC115.03/focused_sequences_msa_trimmed.fasta \
#                                  SPBC115.03/tree/focused_sequences_msa_trimmed.treefile
#
# After building the tree (step 08), these metrics tell you whether the
# alignment + tree actually contain enough signal to be trusted.
#
# How to read the scores (rough rules of thumb):
#   treeness      = internal branch length / total length; higher = more tree-like
#   RCV           = relative composition variability; lower = better
#   treeness/RCV  = the signal-quality ratio; higher = better
#   saturation    = substitution saturation; lower = better (little saturation)
#   pis           = number of parsimony-informative sites (just a count)
#
# Requires: phykit (installed in step 05, into the bioinfo env).

set -e

if ! command -v phykit >/dev/null 2>&1; then
    echo "phykit not found -- install it first (step 05): conda run -n bioinfo pip install phykit"
    exit 1
fi

if [ "$#" -ne 2 ]; then
    echo "usage: bash 09_phylogenetic_qc.sh <trimmed_msa.fasta> <tree.treefile>"
    exit 1
fi

aln=$1
tree=$2
[ -f "$aln" ] || { echo "missing alignment: $aln"; exit 1; }
[ -f "$tree" ] || { echo "missing tree: $tree"; exit 1; }

echo "===== $aln ====="
echo -n "treeness:        "; phykit treeness "$tree"
echo -n "RCV:             "; phykit rcv "$aln"
echo -n "treeness/RCV:    "; phykit treeness_over_rcv -a "$aln" -t "$tree"
echo -n "saturation:      "; phykit saturation -a "$aln" -t "$tree"
echo -n "parsimony sites: "; phykit parsimony_informative_sites "$aln"
