# Step 1: Manual OMA download

Only OMA website downloads are documented here, because scripts are blocked by Cloudflare and the files must be fetched in a browser. All other downloads (tools, PANTHER) are scriptable and go into the .sh file.

## Open the HOG pages

| Gene | HOG page URL |
|---|---|
| SPBC115.03 | https://omabrowser.org/oma/hog/HOG:F0839228.2c.7b/ |
| SPAPB1A10.12c | https://omabrowser.org/oma/hog/HOG:F0849075.5a.15d/ |

## Download these 7 files per level

Use the page's Export/Download menu. For each level, put the files into one folder named `OMA/HOG_xxx_with_N_members/` (N = member count).

| # | Filename | What it is |
|---|---|---|
| 1 | `orthoxml.txt` | orthology structure (OrthoXML) |
| 2 | `speciestree_*.phyloxml` | species tree (no branch lengths) |
| 3 | `fasta.txt` | member sequence list (unaligned; column 3 = sub-HOG path) |
| 4 | `3di-fasta.txt` | 3Di sequences (unaligned) |
| 5 | `protein_sequence_msa.fasta` | official protein MSA (**already aligned**) |
| 6 | `protein_sequence_msa_informed_with_structure.fasta` | structure-informed MSA (**already aligned**) |
| 7 | `structure_alphabet_msa_3di.fasta` | 3Di-alphabet MSA (**already aligned**) |

## Expected folder structure


```
Project_Name/
├── SPBC115.03/
│   └── OMA/
│       ├── HOG_F0839228_with_1102_members/       (7 files)
│       └── HOG_F0839228.2c_with_175_members/     (7 files)
└── SPAPB1A10.12c/
    └── OMA/
        ├── HOG_F0849075_with_617_members/        (7 files)
        ├── HOG_F0849075.5a_with_423_members/     (7 files)
        └── HOG_F0849075.5a.15d_with_226_members/ (7 files)
```

## After download, verify each file

OMA is behind Cloudflare, so a download can occasionally save an HTML challenge page ("Just a moment...") instead of the real data. Check the first line of each file:

| File | Expected first line |
|---|---|
| `orthoxml.txt` | `<?xml ...` or `<orthoXML ...` |
| `speciestree_*.phyloxml` | `<phyloxml ...` |
| `fasta.txt`, `3di-fasta.txt`, `protein_sequence_msa.fasta`, `protein_sequence_msa_informed_with_structure.fasta`, `structure_alphabet_msa_3di.fasta` | starts with `>` |

A file is bad if it starts with `<!DOCTYPE html>` or contains the text "Just a moment" — re-download it.


