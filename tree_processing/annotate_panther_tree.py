#!/usr/bin/env python3
"""
Convert PANTHER NHX annotations into readable node labels.

Takes a Newick tree (already labeled with protein names on the leaves, e.g.
from label_tree.py) whose internal nodes carry NHX annotations like
"[&&NHX:Ev=0>1:S=Ascomycota:ID=AN5]", and replaces each annotation with a
short label:
  S_<species>  for speciation nodes (Ev=0>1)
  D            for duplication nodes (Ev=1>0)
  S / D        for other evidence patterns

The output is a plain Newick tree that Jalview can open.

Usage:
  python3 annotate_panther_tree.py <labeled.nwk> <annotated.nwk>
"""

import re
import sys


def nhx_label(nhx):
    ev = re.search(r'Ev=([01])>([01])', nhx)
    sp = re.search(r'S=([^:]+)', nhx)
    if ev and ev.group(1) == '1' and ev.group(2) == '0':
        return 'D'
    if sp:
        return 'S_' + sp.group(1).strip().replace(' ', '_').replace('-', '_')
    if ev:
        return 'S' if ev.group(2) == '1' else 'D'
    return ''


def sub_with_length(m):
    length, nhx = m.group(1), m.group(2)
    return nhx_label(nhx) + ':' + length


def main():
    if len(sys.argv) != 3:
        print("usage: python3 annotate_panther_tree.py <labeled.nwk> <annotated.nwk>")
        sys.exit(1)

    s = open(sys.argv[1]).read().strip().replace('\n', '').replace('\r', '')

    # internal nodes with a branch length: ":0.556[&&NHX:...]" -> "S_name:0.556"
    s = re.sub(r':([0-9.eE+-]+)\[&&NHX:([^\]]*)\]', sub_with_length, s)
    # root node without a branch length: ")[&&NHX:...]" -> ")S_name"
    s = re.sub(r'\)\[&&NHX:([^\]]*)\]', lambda m: ')' + nhx_label(m.group(1)), s)

    if not s.rstrip().endswith(';'):
        s = s.rstrip() + ';'

    with open(sys.argv[2], 'w') as fh:
        fh.write(s + '\n')

    print(f"wrote {sys.argv[2]} (NHX left: {'[&&NHX' in s})")


if __name__ == "__main__":
    main()
