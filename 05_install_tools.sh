#!/bin/bash
# Step 5: install the analysis tools (alignment + tree + structure).
#   - mafft      : command-line aligner (Homebrew)
#   - Jalview    : GUI alignment viewer (Homebrew cask)
#   - bioinfo env: conda environment with python + clipkit
#     (clipkit pulls biopython/numpy/tqdm as its own dependencies)
#   - iqtree     : maximum-likelihood phylogenetic tree (conda, into bioinfo env; multicore)
#   - phykit     : phylogenetic signal QC (pip, into bioinfo env)
#   - foldseek   : structural alignment / TM-scores (conda, into bioinfo env)
#   - ChimeraX   : 3D structure viewer (manual download, UCSF site)
# Run section by section. Skip a section if the tool is already installed.

set -e

# ---------------------------------------------------------------
# Section 1. mafft (sequence alignment)
# ---------------------------------------------------------------
# Check if it is already installed:
if command -v mafft >/dev/null 2>&1; then
    echo "mafft already installed: $(command -v mafft)"
else
    brew install mafft
fi

# ---------------------------------------------------------------
# Section 2. Jalview (GUI viewer; also usable headless from CLI)
# ---------------------------------------------------------------
# Check if the app is already installed:
if [ -d "/Applications/Jalview.app" ]; then
    echo "Jalview already installed"
else
    brew install --cask jalview
fi

# ---------------------------------------------------------------
# Section 3. bioinfo conda environment (python + clipkit)
# ---------------------------------------------------------------
# CONDA_SUBDIR=osx-arm64 forces native Apple Silicon packages.
if conda env list | grep -qw bioinfo; then
    echo "bioinfo env already exists"
else
    CONDA_SUBDIR=osx-arm64 conda create -n bioinfo python=3.12 clipkit \
        -c conda-forge -c bioconda -y --override-channels
fi

# Activate it in every new terminal before analysis:
#   conda activate bioinfo
# Verify (should print "arm64"):
#   python3 -c "import platform; print(platform.machine())"

# ---------------------------------------------------------------
# Section 4. iqtree (maximum-likelihood phylogenetic tree, multicore)
# ---------------------------------------------------------------
# Install into the bioinfo env (created in Section 3).
# The conda build is multicore (OpenMP); the Homebrew build is single-core,
# which errors on "-nt AUTO". Requires the bioinfo env to be active to use.
CONDA_SUBDIR=osx-arm64 conda install -n bioinfo iqtree -c bioconda -c conda-forge -y --override-channels

# ---------------------------------------------------------------
# Section 5. phykit (phylogenetic signal QC)
# ---------------------------------------------------------------
# Install into the bioinfo env (created in Section 3).
conda run -n bioinfo pip install phykit

# ---------------------------------------------------------------
# Section 6. foldseek (structural alignment / TM-scores)
# ---------------------------------------------------------------
# Install into the bioinfo env (created in Section 3).
# --override-channels skips the configured conda mirrors (which can time out);
# CONDA_SUBDIR forces the native Apple Silicon build.
CONDA_SUBDIR=osx-arm64 conda install -n bioinfo foldseek -c conda-forge -c bioconda -y --override-channels

# ---------------------------------------------------------------
# Section 7. ChimeraX (3D structure viewer; no Homebrew cask)
# ---------------------------------------------------------------
# Manual download from the UCSF site (not available via Homebrew):
#   https://www.cgl.ucsf.edu/chimerax/download.html
# Download the macOS .dmg, open it, and drag ChimeraX.app into /Applications.
if ls /Applications/ChimeraX*.app >/dev/null 2>&1; then
    echo "ChimeraX already installed"
else
    echo "ChimeraX not found -- download from https://www.cgl.ucsf.edu/chimerax/download.html"
fi
