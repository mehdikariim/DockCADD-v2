# DockCADD-v3

DockCADD-v3 is an end-to-end docking workflow built around RDKit, P2Rank, AutoDock Vina, Open Babel, and optional PDBFixer cleanup.

It supports:

- SMILES input
- ligand CSV files
- ligand SDF files
- ligand PDB files
- downloaded receptors from PDB IDs
- local receptor PDB files
- redocking from a known ligand structure
- matrix runs with CSV and Excel reporting

The installable package name is `dockcadd`.

## Quick Start

Single target from SMILES:

```bash
dockcadd dock --pdb-id 5TZ1 --smiles "CCO" "CCN" --output-root outputs/5TZ1
```

Single target from an SDF file:

```bash
dockcadd dock --pdb-id 5TZ1 --ligands-sdf ligands.sdf --output-root outputs/5TZ1
```

Single target from local receptor and ligand PDB files:

```bash
dockcadd dock --receptor-pdb receptors/target.pdb --ligands-pdb ligands/ligand1.pdb ligands/ligand2.pdb --output-root outputs/local_target
```

Redocking:

```bash
dockcadd redock --pdb-id 5TZ1 --ligand-file ligands/reference.sdf --output-root outputs/redock_5TZ1
```

Matrix run:

```bash
dockcadd matrix --config-json examples/chebbak_matrix.json
```

## Inputs

### Ligands

DockCADD-v3 accepts ligands as:

- `--smiles`
- `--ligands-csv` with a SMILES column
- `--ligands-sdf`
- `--ligands-pdb`

### Receptors

DockCADD-v3 accepts receptors as:

- `--pdb-id`
- `--receptor-pdb`

Exactly one receptor source must be provided for `dock` and `redock`.

## Installation

### Linux

Clone and enter the repository:

```bash
git clone https://github.com/mehdikariim/DockCADD-v2.git
cd DockCADD-v2
```

Install system dependencies and the package:

```bash
bash scripts/setup.sh
pip install -e .
```

Run:

```bash
dockcadd dock --pdb-id 5TZ1 --smiles "CCO" --output-root outputs/5TZ1
```

### Windows With WSL

Clone the repository in Windows or in WSL, then open WSL in the repository directory:

```bash
git clone https://github.com/mehdikariim/DockCADD-v2.git
cd DockCADD-v2
bash scripts/setup.sh
pip install -e .
```

Run:

```bash
dockcadd dock --pdb-id 5TZ1 --smiles "CCO" --output-root outputs/5TZ1
```

### Native Windows

Requirements:

- Python 3.10+
- Java
- Open Babel
- AutoDock Vina

Clone and bootstrap:

```powershell
git clone https://github.com/mehdikariim/DockCADD-v2.git
cd DockCADD-v2
powershell -ExecutionPolicy Bypass -File scripts/setup_windows.ps1
pip install -e .
```

Run:

```powershell
dockcadd dock --pdb-id 5TZ1 --smiles "CCO" --output-root outputs\5TZ1
```

### Google Colab

Setup:

```python
!git clone https://github.com/mehdikariim/DockCADD-v2.git
%cd DockCADD-v2
!bash scripts/setup.sh
!pip install -e .
```

Run:

```python
!dockcadd dock --pdb-id 5TZ1 --smiles "CCO" --output-root outputs/5TZ1
```

## Commands

### 1. Dock from SMILES

```bash
dockcadd dock --pdb-id 5TZ1 --smiles "CCO" "CCN" --output-root outputs/5TZ1
```

### 2. Dock from a CSV file

```bash
dockcadd dock --pdb-id 5TZ1 --ligands-csv ligands.csv --smiles-column SMILES --output-root outputs/5TZ1
```

### 3. Dock from one or more SDF files

```bash
dockcadd dock --pdb-id 5TZ1 --ligands-sdf ligands.sdf more_ligands.sdf --output-root outputs/5TZ1
```

### 4. Dock from one or more ligand PDB files

```bash
dockcadd dock --pdb-id 5TZ1 --ligands-pdb ligands/a.pdb ligands/b.pdb --output-root outputs/5TZ1
```

### 5. Dock against a local receptor PDB

```bash
dockcadd dock --receptor-pdb receptors/target.pdb --ligands-csv ligands.csv --output-root outputs/local_target
```

### 6. Redocking

```bash
dockcadd redock --pdb-id 5TZ1 --ligand-file ligands/reference.sdf --output-root outputs/redock_5TZ1
```

For a local receptor:

```bash
dockcadd redock --receptor-pdb receptors/target.pdb --ligand-file ligands/reference.pdb --output-root outputs/redock_local
```

### 7. Matrix Mode

```bash
dockcadd matrix --config-json examples/chebbak_matrix.json
```

## Matrix Config

See:

- `examples/chebbak_matrix.json`

Paths in the config can be relative. They are resolved relative to the config file location.

Each target can define either:

- `"pdb_id": "5TZ1"`
- `"receptor_pdb": "receptors/target.pdb"`

## Output

Each target output folder contains:

- receptor `.pdb`
- receptor `.pdbqt`
- ligand `.pdb`
- ligand `.pdbqt`
- Vina logs
- docked poses
- `docking_results.txt`

Matrix mode also writes:

- Excel workbook
- CSV summary
- `DockingResults` sheet
- `TargetSummary` sheet
- `CompoundResolution` sheet

## Notes

- If `pdbfixer` is installed, DockCADD-v3 uses it to clean receptor structures before docking.
- If `pdbfixer` is not installed, DockCADD-v3 falls back to conservative PDB cleanup.
- On Linux or Colab, `scripts/setup.sh` installs the core stack and attempts `openmm`; `pdbfixer` may still require a conda-based install depending on the environment.
- On Windows, full `pdbfixer` support is best handled through WSL or a conda environment.
- `redock` preserves the supplied ligand structure file and docks that known ligand back into the receptor workflow.

## Citation

If you use DockCADD-v3 in scientific work, cite:

- Karim El Mehdi et al., *Scientific African* (2025). DOI: https://doi.org/10.1016/j.sciaf.2025.e02581
