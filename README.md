# Unstudied Genes — Crick 2026

A reproducible, from-scratch bioinformatics workflow characterising two
*Schizosaccharomyces pombe* "priority unstudied genes".

## Target genes

| Gene | UniProt | Putative function |
|---|---|---|
| SPBC115.03 | O42896 | gfo/idh/mocA oxidoreductase family — NAD(P)-dependent sugar dehydrogenase |
| SPAPB1A10.12c | Q9HDX8 | D-arabinono-1,4-lactone oxidase (ALO) — FAD-dependent, erythroascorbate biosynthesis |

## Workflow (steps 01–14)

| Step | Script | What it does |
|---|---|---|
| 01 | `01_manual_oma_download.md` | Download OMA HOG / sub-HOG |
| 02 | `02_panther_download.sh` | Download PANTHER family |
| 03 | `03_screen_focused_genes.py` | Build the focused gene set |
| 04 | `04_fetch_sequences.sh` | Fetch protein sequences |
| 05 | `05_install_tools.sh` | Install all tools |
| 06 | `06_pairwise_percent_id.py` | Pairwise sequence identity |
| 07 | `07_align_and_trim.sh` | mafft alignment + ClipKIT trimming |
| 08 | `08_build_tree.sh` | IQ-TREE phylogeny |
| 09 | `09_phylogenetic_qc.sh` | phykit tree QC |
| 10 | `10_conservation.py` | Conservation scores (JSD) |
| 11 | `11_download_alphafold.sh` | Download AlphaFold structures |
| 12 | `12_structural_alignment.sh` | Foldseek TM-scores |
| 13 | `13_candidate_residues.py` | Candidate-residue list (conservation + annotation + literature) |
| 14 | `14_chimerax_script.py` | ChimeraX colouring script |

## Repository layout

- `<GENE>/` — per-gene results: candidate residues, annotations, TM-scores, MSA, etc.
- `tree_processing/` — tree-conversion utilities

Large regenerable data (OMA/PANTHER MSAs, AlphaFold PDBs, IQ-TREE outputs) is
**not** committed.