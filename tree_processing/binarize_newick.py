#!/usr/bin/env python3
"""
Resolve multifurcations in a Newick tree into bifurcations.

Jalview has a rendering bug with >=4-way polytomies, so trees are binarized:
each polytomy is split by adding dummy nodes with a tiny 0.0001 branch length.

Usage:
  python3 binarize_newick.py <in.nwk> <out.nwk>
"""

import sys


def tokenize(s):
    toks = []
    i = 0
    while i < len(s):
        c = s[i]
        if c in '(),':
            toks.append(c)
            i += 1
        elif c == ':':
            j = i + 1
            while j < len(s) and s[j] not in '(),':
                j += 1
            toks.append(s[i:j])
            i = j
        elif c == '[':
            j = s.index(']', i)
            toks.append(s[i:j + 1])
            i = j + 1
        else:
            j = i + 1
            while j < len(s) and s[j] not in '(),:[':
                j += 1
            toks.append(s[i:j])
            i = j
    return toks


class Node:
    def __init__(self):
        self.label = None
        self.bl = None
        self.children = []
        self.nhx = None
        self.dummy = False


def parse(toks, pos):
    n = Node()
    if toks[pos] == '(':
        pos += 1
        while toks[pos] != ')':
            c, pos = parse(toks, pos)
            n.children.append(c)
            if toks[pos] == ',':
                pos += 1
        pos += 1
    if pos < len(toks) and toks[pos] not in (':', ',', ')') \
            and not toks[pos].startswith('[') and not toks[pos].startswith(':'):
        n.label = toks[pos]
        pos += 1
    if pos < len(toks) and toks[pos].startswith(':'):
        n.bl = toks[pos]
        pos += 1
    if pos < len(toks) and toks[pos].startswith('['):
        n.nhx = toks[pos]
        pos += 1
    return n, pos


def binarize(n):
    kids = list(n.children)
    while len(kids) > 2:
        d = Node()
        d.dummy = True
        d.children = [kids[0], kids[1]]
        kids = [d] + kids[2:]
    n.children = kids
    for c in n.children:
        binarize(c)


def serialize(n):
    if n.dummy:
        return '(' + ','.join(serialize(c) for c in n.children) + '):0.0001'
    if n.children:
        inner = ','.join(serialize(c) for c in n.children)
        return '(' + inner + ')' + (n.label or '') + (n.bl or '') + (n.nhx or '')
    return (n.label or '') + (n.bl or '') + (n.nhx or '')


def main():
    if len(sys.argv) != 3:
        print("usage: python3 binarize_newick.py <in.nwk> <out.nwk>")
        sys.exit(1)
    s = open(sys.argv[1]).read().replace('\n', '').replace('\r', '').strip()
    if s.endswith(';'):
        s = s[:-1]
    root, _ = parse(tokenize(s), 0)
    binarize(root)
    out = serialize(root) + ';\n'
    open(sys.argv[2], 'w').write(out)
    print(f"wrote {sys.argv[2]}")


if __name__ == "__main__":
    main()
