#!/usr/bin/env python3
"""
Label the leaves of a PANTHER tree.

PANTHER's tree.tree is a two-part file: line 1 is a Newick tree whose leaves
are internal node IDs ("AN6"), and the remaining lines are "AN6:<label>"
mappings to the protein labels. This script replaces each leaf ID with its
label, producing a Newick tree that Jalview can read.

Usage:
  python3 label_tree.py <tree.tree> <out.nwk>
"""

import re
import sys


def main():
    if len(sys.argv) != 3:
        print("usage: python3 label_tree.py <tree.tree> <out.nwk>")
        sys.exit(1)

    lines = open(sys.argv[1]).read().splitlines()
    newick = lines[0]

    mapping = {}
    for line in lines[1:]:
        m = re.match(r'^(AN\d+):(.+?);?\s*$', line)
        if m:
            mapping[m.group(1)] = m.group(2)

    def replace_leaf(m):
        an = m.group(1)
        return mapping[an] + ':' if an in mapping else m.group(0)

    labeled = re.sub(r'(?<![A-Za-z0-9])(AN\d+):(?=[0-9.]+)', replace_leaf, newick)

    with open(sys.argv[2], 'w') as fh:
        fh.write(labeled + '\n')

    remaining = re.findall(r'(?<![A-Za-z0-9])(AN\d+):(?=[0-9.]+)', labeled)
    print(f"mapped {len(mapping)} leaves, {len(remaining)} unresolved -> {sys.argv[2]}")


if __name__ == "__main__":
    main()
