#!/usr/bin/env python3
"""
Convert an OMA species tree (speciestree_*.phyloxml) into Newick.

Jalview cannot read phyloXML, so this rebuilds the species tree as a Newick
string, using the scientific names as node labels and binarizing any polytomy.

Usage:
  python3 phyloxml2newick.py <OMA_folder>
    <OMA_folder>  directory containing speciestree_*.phyloxml

Output:
  <OMA_folder>/OMA_<name>_species_tree.nwk
"""

import glob
import re
import sys
import xml.etree.ElementTree as ET

PHYLO_NS = {'p': 'http://www.phyloxml.org'}


def clade_newick(el):
    kids = el.findall('p:clade', PHYLO_NS)

    label = ''
    tax = el.find('p:taxonomy/p:scientific_name', PHYLO_NS)
    if tax is not None and tax.text:
        label = tax.text.replace(' ', '_')
        for ch in '()/:,;\'"[]':
            label = label.replace(ch, '_')

    bl = el.find('p:branch_length', PHYLO_NS)
    suffix = ':' + bl.text if (bl is not None and bl.text) else ':1'

    if kids:
        if len(kids) == 1:
            return clade_newick(kids[0])
        parts = [clade_newick(k) for k in kids]
        while len(parts) > 2:
            parts = [f'({parts[0]},{parts[1]}):0.0001'] + parts[2:]
        return '(' + ','.join(parts) + ')' + label + suffix
    return (label or 'leaf') + suffix


def main():
    if len(sys.argv) != 2:
        print("usage: python3 phyloxml2newick.py <OMA_folder>")
        sys.exit(1)
    base = sys.argv[1]

    path = glob.glob(f'{base}/speciestree_*.phyloxml')[0]
    root = ET.parse(path).getroot()

    phyl = root.find('p:phylogeny', PHYLO_NS)
    top = phyl.findall('p:clade', PHYLO_NS)
    parts = [clade_newick(c) for c in top]
    if len(parts) == 1:
        nwk = parts[0] + ';'
    else:
        while len(parts) > 2:
            parts = [f'({parts[0]},{parts[1]}):0.0001'] + parts[2:]
        nwk = '(' + ','.join(parts) + ')root:1;'

    m = re.search(r'speciestree_(\w+)\.phyloxml', path)
    tree_name = m.group(1) if m else 'species'
    out = f'{base}/OMA_{tree_name}_species_tree.nwk'
    with open(out, 'w') as fh:
        fh.write(nwk + '\n')
    print(f'species tree written: {out}')


if __name__ == "__main__":
    main()
