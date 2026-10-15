#!/bin/bash
# Step 4: download sequences for the focused genes from UniProt.
#
# One input file per gene: <GENE>/focused_genes.txt
#         One UniProt accession per line. An optional label may follow the accession, separated by whitespace:
#             O42896   SPBC115.03_pombe
#             Q9UT60
#         Blank lines and lines starting with "#" are ignored.
# One output file per gene: <GENE>/focused_sequences.fasta  (unaligned)

set -e

# Download the focused sequences for one gene.
fetch_gene () {
    local gene=$1
    local input="${gene}/focused_genes.txt"
    local output="${gene}/focused_sequences.fasta"

    if [ ! -f "$input" ]; then
        echo "missing input file: $input"
        exit 1
    fi

    : > "$output"

    while read -r accession label; do
        # Skip blank lines and comment lines.
        [ -z "$accession" ] && continue
        case "$accession" in
            \#*) continue ;;
        esac

        # Build the FASTA header: ">label|accession" if a label is given, otherwise just ">accession".
        if [ -n "$label" ]; then
            header=">${label}|${accession}"
        else
            header=">${accession}"
        fi

        echo "$header" >> "$output"

        # UniProt returns the sequence in FASTA format; we drop its own header line and keep only the sequence.
        curl -sL --max-time 60 "https://rest.uniprot.org/uniprotkb/${accession}.fasta" \
            | grep -v '^>' >> "$output"
    done < "$input"

    echo "wrote $output"
}

fetch_gene SPBC115.03
fetch_gene SPAPB1A10.12c
