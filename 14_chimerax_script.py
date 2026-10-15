#!/usr/bin/env python3
"""
Step 14: generate a ChimeraX command script (.cxc) that colours the candidate
residues (from step 13's candidate_residues.csv) by source, so the 3D mapping
is reproducible.

Colours:
  annotation   red     (known UniProt residue, e.g. the FAD-histidine)
  literature   orange  (from papers, mapped from a homolog)
  conservation yellow  (predicted from the MSA)

Usage:
  python3 14_chimerax_script.py <GENE>

Output:
  <GENE>/candidate_residues.cxc

Then, with ChimeraX installed, run it in the GUI (headless --nogui / --offscreen
cannot render/save images on macOS because OpenGL is unavailable):
  open <GENE>/candidate_residues.cxc
(The trailing "save ... .png" line works in the GUI, not headless.)
"""

import argparse
import csv
import os

GOI_ACCESSION = {"SPBC115.03": "O42896", "SPAPB1A10.12c": "Q9HDX8"}
COLORS = {"annotation": "red", "literature": "orange", "conservation": "yellow"}

# Absolute paths, so the .cxc works no matter where ChimeraX is launched from
# (macOS GUI apps do not inherit the shell's working directory).
BASE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gene", choices=["SPBC115.03", "SPAPB1A10.12c"])
    args = ap.parse_args()
    gene = args.gene
    goi = GOI_ACCESSION[gene]

    # Group residue numbers by source.
    by_source = {"annotation": [], "literature": [], "conservation": []}
    for r in csv.reader(open(f"{gene}/candidate_residues.csv")):
        if r and r[0] != "goi_residue":
            by_source[r[2]].append(r[0])

    pdb = os.path.join(BASE, gene, "alphafold", f"{goi}.pdb")
    png = os.path.join(BASE, gene, f"{goi}_candidates.png")

    lines = [f"open {pdb}"]
    for source, residues in by_source.items():
        if residues:
            lines.append(f"color :{','.join(residues)} {COLORS[source]}")
    lines.append(f"save {png}")

    out = f"{gene}/candidate_residues.cxc"
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"wrote {out}")
    for line in lines:
        print(f"  {line}")


if __name__ == "__main__":
    main()
