# Tree processing utilities

These scripts convert the raw OMA and PANTHER tree files (which Jalview cannot
open directly) into plain Newick trees that Jalview can read. They were written
in the earlier Crick2026 work and copied here **verbatim** so the conversion is
reproducible.

## Why the raw trees do not open in Jalview

- **PANTHER `tree.tree`** is Newick with NHX annotations (`[&&NHX:...]`), and
  its leaf labels are internal node IDs (`AN6`, `AN8`, ...), not protein names.
- **OMA `orthoxml.txt`** is OrthoXML (XML), which Jalview cannot parse.
- Jalview has a rendering bug with >=4-way polytomies, so trees must be
  binarized (multifurcations resolved to bifurcations).

## The conversion pipeline

### PANTHER (`tree.tree` -> Newick)

1. `label_tree.py` — replace the `AN#` leaf IDs with protein labels.
   Input: the PANTHER `tree.tree` file **directly** — its line 1 is the Newick
   and its later lines already carry the `AN#:<label>` mappings, so no separate
   mapping file is needed.
2. `annotate_panther_tree.py` — turn the NHX annotations into readable labels:
   `S_<species>` for speciation nodes, `D` for duplication nodes.
   Output: `*_tree_annotated.nwk` (Jalview-openable).

### OMA (`orthoxml.txt` -> Newick)

- `orthoxml2newick_v2.py` — parse the OrthoXML gene tree (the
  `orthologGroup` / `paralogGroup` structure), map gene ids -> protIds ->
  UniProt accessions (via `fasta.txt`), and write
  `OMA_*_official_structure_tree.nwk`. It also binarizes the tree. Supports
  pruning to a member list or selecting one sub-HOG.

### Species tree (`speciestree_*.phyloxml` -> Newick)

- `phyloxml2newick.py` — convert the OMA species tree (phyloXML) to Newick.

### Binarization (fixes the Jalview polytomy bug)

- `binarize_newick.py` — take any Newick, resolve >=3-way multifurcations into
  bifurcations using dummy nodes with a tiny `:0.0001` branch length.

## Usage

Each script takes input/output paths as command-line arguments. A few of them
still carry the hardcoded default paths from the original Crick2026 run
(e.g. `annotate_panther_tree.py`, `phyloxml2newick.py`); pass explicit
arguments or edit those defaults before re-running.

## Where these fit in the Crick2026rep workflow

Steps 01 (OMA) and 02 (PANTHER) download the raw data. Run these scripts on the
downloaded trees whenever you want to view them in Jalview. They are kept here
as standalone utilities rather than numbered steps, because they are not part
of the focused-gene analysis pipeline (06-13).
