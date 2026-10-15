#!/bin/bash
# Step 2: download PANTHER family data.
# The PANTHER data server has no Cloudflare, so curl works directly.
# After running, check that each file is a real data file (not an HTML error page).

set -e

PANTHER_BASE=https://data.pantherdb.org/PANTHER19.0/books

# Download the 5 family files into <gene>/PANTHER/<family_id>/.
download_panther () {
    local gene=$1
    local fam=$2
    mkdir -p "${gene}/PANTHER/${fam}"
    for f in tree.tree cluster.fasta attr.tab cluster.ortholog tree.mia; do
        curl -L -o "${gene}/PANTHER/${fam}/${f}" "${PANTHER_BASE}/${fam}/${f}"
    done
}

# download_panther SPBC115.03 PTHR43708   # SPBC115.03 family
# download_panther SPAPB1A10.12c PTHR43762   # SPAPB1A10.12c family
# download_panther SPBC115.03 PTHR42840   # YMR315W
# download_panther SPBC115.03 PTHR43818   # HUMAN GFOD1 and GFOD2
# download_panther SPBC115.03 PTHR22604   # HUMAN DHDH
# download_panther SPAPB1A10.12c PTHR13878   # DprE1
