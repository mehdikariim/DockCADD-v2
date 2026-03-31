# DockCADD-v3

DockCADD-v3 is an end-to-end molecular docking workflow built for practical screening runs on Linux, Windows/WSL, and Google Colab. It automates ligand preparation, receptor download and cleanup, pocket detection with P2Rank, docking with AutoDock Vina, and tabular result export.

## What changed in this refactor

- unified the codebase under a single package name: `dockcadd`
- kept compatibility wrappers for older `src.cadock` and `caddock.docking` imports
- removed PyMOL as a required runtime dependency
- exposed a CLI entrypoint: `dockcadd`
- improved ligand preparation robustness with an Open Babel fallback
- made output handling safer and easier to reuse across multiple targets

## Core workflow

1. Accept ligands as SMILES from the CLI or a CSV file
2. Download the receptor structure from the PDB
3. Remove `HETATM` records from the receptor
4. Predict the top binding pocket with P2Rank
5. Build the docking box from the predicted pocket residues
6. Convert receptor and ligands to PDBQT with Open Babel
7. Dock ligands with AutoDock Vina
8. Write ranked docking scores and docking poses to disk

## Requirements

- Python 3.10+
- Open Babel available in `PATH`
- AutoDock Vina available in `PATH`
- Java available for P2Rank
- an extracted `p2rank_2.4.2` directory in the repository root, or pass `--p2rank-dir`

## Installation

### Python package

```bash
pip install -r requirements.txt
pip install -e .
```

### Linux / Google Colab bootstrap

```bash
bash scripts/setup.sh
```

### Windows

Recommended:

- run DockCADD inside WSL, or
- install Python, Java, Open Babel, and AutoDock Vina natively and keep them in `PATH`
- then bootstrap:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/setup_windows.ps1
```

## CLI examples

Single target from direct SMILES:

```bash
dockcadd dock --pdb-id 5TZ1 --smiles "CCO" "CCN" --output-root results/5TZ1
```

From a CSV file:

```bash
dockcadd dock --pdb-id 5TZ1 --ligands-csv ligands.csv --smiles-column SMILES --output-root results/5TZ1
```

Matrix run from JSON config:

```bash
dockcadd matrix --config-json examples/chebbak_matrix.json
```

## Python example

```python
from pathlib import Path

from dockcadd import perform_docking, run_matrix_from_config

perform_docking(
    ["CCO", "CCN"],
    "5TZ1",
    output_root="results/5TZ1",
)

run_matrix_from_config(Path("examples/chebbak_matrix.json"))
```

## Outputs

Each target run writes:

- cleaned receptor files
- ligand PDB/PDBQT files
- Vina log files
- docked pose files
- `docking_results.txt`

Matrix runs also write:

- `DockingResults` sheet
- `TargetSummary` sheet
- `CompoundResolution` sheet
- CSV export of the full matrix

## Example dataset

- `examples/chebbak_matrix.json` reproduces the Artemisia docking matrix used in this workspace.

## Citation

If you use DockCADD-v3 in scientific work, cite:

- Karim El Mehdi et al., *Scientific African* (2025), DOI: https://doi.org/10.1016/j.sciaf.2025.e02581

## Current scope

DockCADD is intended to be a practical automated docking workflow. It is not a replacement for full receptor curation, protonation-state studies, induced-fit modeling, or post-docking rescoring pipelines.
