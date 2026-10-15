#!/bin/bash
# Step 11: download AlphaFold predicted structures for the focused genes.
#
# Input : <GENE>/focused_sequences.fasta   (one UniProt accession per header,
#                                           from step 04)
# Output: <GENE>/alphafold/<ACCESSION>.pdb  (one AlphaFold model per accession)
#
# Uses the AlphaFold DB API to find the current PDB URL, then downloads it.
# Requires curl (already on macOS).

set -e

download_one () {
    local accession=$1
    local outdir=$2
    local out="${outdir}/${accession}.pdb"

    [ -f "$out" ] && return  # already downloaded

    # Ask AlphaFold DB where this accession's model lives.
    local pdb_url
    pdb_url=$(curl -sL --max-time 60 "https://alphafold.ebi.ac.uk/api/prediction/${accession}" \
        | python3 -c 'import sys, json; d = json.load(sys.stdin); print(d[0]["pdbUrl"] if d else "")')

    if [ -z "$pdb_url" ]; then
        echo "no AlphaFold model for $accession"
        return
    fi

    # The AlphaFold DB returns 403 to curl's default user agent, so send a
    # browser-like User-Agent. Then check the file is a real PDB (not an
    # error page) before keeping it.
    curl -sL -A "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36" \
        --max-time 120 -o "$out" "$pdb_url"

    if grep -q '^HEADER\|^ATOM' "$out"; then
        echo "downloaded $accession -> $out"
    else
        echo "download failed for $accession (not a PDB); removing $out"
        rm -f "$out"
    fi
}

for gene in SPBC115.03 SPAPB1A10.12c; do
    fasta="${gene}/focused_sequences.fasta"
    outdir="${gene}/alphafold"
    [ -f "$fasta" ] || { echo "missing input: $fasta"; continue; }
    mkdir -p "$outdir"

    # Extract the UniProt accession from each header (the part after "|").
    grep '^>' "$fasta" | sed 's/^>//' | while read -r header; do
        accession="${header##*|}"
        download_one "$accession" "$outdir"
    done
done
