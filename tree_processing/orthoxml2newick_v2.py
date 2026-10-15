#!/usr/bin/env python3
"""
Convert an OMA HOG orthoxml.txt tree into Newick.

OMA stores its gene tree inside orthoxml.txt as nested <orthologGroup> /
<paralogGroup> elements (speciation vs. duplication nodes) ending in <geneRef>
leaves. Jalview cannot read OrthoXML, so this script rebuilds the tree as a
Newick string: leaves are labeled "<protId>|<UniProtAC>" (from fasta.txt),
internal nodes are labeled "S_<species>" / "D" (from the species tree), and the
tree is binarized to avoid Jalview's polytomy rendering bug.

Usage:
  python3 orthoxml2newick_v2.py <OMA_folder> [--subhog <id>] [--members <fasta>] [--label <name>]
    <OMA_folder>  directory containing orthoxml.txt, fasta.txt, speciestree_*.phyloxml
    --subhog      only output the sub-HOG whose id starts with "HOG:<id>_"
    --members     prune to the sequences whose headers appear in this fasta
    --label       prefix for the output file name (default: the subhog id or "HOG")

Output (two per run):
  <OMA_folder>/OMA_<label>_official_structure_tree.nwk           (leaves = protId|AC)
  <OMA_folder>/OMA_<label>_official_structure_tree_msa_sync.nwk  (leaves = protId)
"""

import argparse
import glob
import sys
import xml.etree.ElementTree as ET

ORTHO_NS = {'o': 'http://orthoXML.org/2011/'}
PHYLO_NS = {'p': 'http://www.phyloxml.org'}


class Node:
    def __init__(self, label=None, kind='leaf'):
        self.label = label
        self.kind = kind  # 'leaf', 'S' (speciation), 'D' (duplication), 'dummy'
        self.children = []


def read_fasta_headers(fasta_path):
    """Return {protId: uniprot_ac} from an OMA fasta.txt."""
    prot2ac = {}
    for line in open(fasta_path):
        if line.startswith('>'):
            parts = [p.strip() for p in line[1:].strip().split('|')]
            if len(parts) >= 2:
                prot2ac[parts[0]] = parts[1]
    return prot2ac


def build_taxonomy_map(oma_folder):
    """Map taxon id -> species name from the species tree phyloXML."""
    path = glob.glob(f'{oma_folder}/speciestree_*.phyloxml')[0]
    root = ET.parse(path).getroot()
    tax2name = {}
    for el in root.iter():
        if el.tag.split('}')[-1] == 'taxonomy':
            tid = el.find('p:id', PHYLO_NS)
            sn = el.find('p:scientific_name', PHYLO_NS)
            if tid is not None and sn is not None:
                tax2name[tid.text] = sn.text
    return tax2name


def build(el, gene2prot, tax2name):
    """Recursively rebuild the tree from an orthoXML element."""
    tag = el.tag.split('}')[-1]
    if tag == 'geneRef':
        return Node(gene2prot.get(el.get('id'), el.get('id')), 'leaf')

    n = Node(None, 'S' if tag == 'orthologGroup' else 'D')
    taxid = el.get('taxonId', '')
    if taxid and taxid in tax2name:
        n.label = tax2name[taxid]
    for c in el:
        if c.tag.split('}')[-1] in ('orthologGroup', 'paralogGroup', 'geneRef'):
            n.children.append(build(c, gene2prot, tax2name))
    return n


def prune(n, keep):
    """Keep only the leaves whose protId is in `keep`; collapse single-child nodes."""
    if n.kind == 'leaf':
        return n if n.label in keep else None
    kids = [prune(c, keep) for c in n.children]
    kids = [k for k in kids if k is not None]
    if not kids:
        return None
    if len(kids) == 1:
        return kids[0]
    n.children = kids
    return n


def binarize(n):
    """Resolve polytomies by inserting dummy nodes (fixes the Jalview bug)."""
    kids = list(n.children)
    while len(kids) > 2:
        dummy = Node(None, 'dummy')
        dummy.children = [kids[0], kids[1]]
        kids = [dummy] + kids[2:]
    n.children = kids
    for c in n.children:
        binarize(c)


def serialize(n, leaf_style, prot2ac):
    if n.kind == 'leaf':
        lab = n.label
        if leaf_style == 'full' and lab in prot2ac:
            lab = f"{lab}|{prot2ac[lab].replace(':', '_')}"
        return lab + ':1'
    if n.kind == 'dummy':
        inner = ','.join(serialize(c, leaf_style, prot2ac) for c in n.children)
        return f'({inner}):0.0001'
    inner = ','.join(serialize(c, leaf_style, prot2ac) for c in n.children)
    if n.label:
        tname = n.label.replace(' ', '_').replace(':', '_')
        lab = f'{n.kind}_{tname}'
    else:
        lab = n.kind
    return f'({inner}){lab}:1'


def main():
    parser = argparse.ArgumentParser(description="Convert OMA orthoxml tree to Newick.")
    parser.add_argument("folder", help="OMA HOG folder with orthoxml.txt, fasta.txt, speciestree_*.phyloxml")
    parser.add_argument("--subhog", help="only output this sub-HOG (id starts with HOG:<subhog>_)")
    parser.add_argument("--members", help="fasta file whose headers define the leaves to keep")
    parser.add_argument("--label", help="output file name prefix")
    args = parser.parse_args()

    root = ET.parse(f'{args.folder}/orthoxml.txt').getroot()

    gene2prot = {}
    for sp in root.findall('o:species', ORTHO_NS):
        for db in sp.findall('o:database', ORTHO_NS):
            for genes in db.findall('o:genes', ORTHO_NS):
                for g in genes.findall('o:gene', ORTHO_NS):
                    gene2prot[g.get('id')] = g.get('protId')

    try:
        prot2ac = read_fasta_headers(f'{args.folder}/fasta.txt')
    except FileNotFoundError:
        prot2ac = {}

    tax2name = build_taxonomy_map(args.folder)

    groups_el = root.find('o:groups', ORTHO_NS)
    top = [c for c in groups_el if c.tag.split('}')[-1] in ('orthologGroup', 'paralogGroup')]
    if len(top) != 1:
        sys.exit(f'expected 1 top-level group, found {len(top)}')

    t = build(top[0], gene2prot, tax2name)
    label = args.label or 'HOG'

    if args.subhog:
        prefix = f'HOG:{args.subhog}_'
        found = None
        for el in root.iter():
            tag = el.tag.split('}')[-1]
            if tag in ('orthologGroup', 'paralogGroup') and (el.get('id') or '').startswith(prefix):
                found = el
                break
        if found is None:
            sys.exit(f'sub-HOG not found: {args.subhog}')
        t = build(found, gene2prot, tax2name)
        label = args.subhog
    elif args.members:
        keep = set()
        for line in open(args.members):
            line = line.strip()
            if line.startswith('>'):
                keep.add(line[1:].split('|')[0].strip())
        t = prune(t, keep)
        if t is None:
            sys.exit('pruning removed every leaf')

    binarize(t)

    for style, suffix in [('full', 'official_structure_tree'),
                          ('plain', 'official_structure_tree_msa_sync')]:
        nwk = serialize(t, style, prot2ac)
        outname = f'{args.folder}/OMA_{label}_{suffix}.nwk'
        with open(outname, 'w') as fh:
            fh.write(nwk + ';\n')
        print(f'wrote {outname}')


if __name__ == "__main__":
    main()
