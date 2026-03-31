from __future__ import annotations

from pathlib import Path
import os
import shutil
import subprocess
import sys
import tempfile

from Bio.PDB import PDBList, PDBParser
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from .config import DockingRunConfig


def generate_minimized_pdb(smiles: str, pdb_filename: str | Path, *, obabel: str = "obabel") -> bool:
    """Generate a 3D ligand structure and write it as PDB."""
    output_path = Path(pdb_filename)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES string: {smiles}")
        return False

    mol = Chem.AddHs(mol)
    try:
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        embed_status = AllChem.EmbedMolecule(mol, params)
        if embed_status != 0:
            embed_status = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
        if embed_status != 0:
            raise ValueError("RDKit embedding failed")
    except Exception:
        mol = None

    try:
        if mol is not None:
            if AllChem.MMFFHasAllMoleculeParams(mol):
                mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
                if AllChem.MMFFOptimizeMolecule(mol, mmff_props, maxIters=500) != 0:
                    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            else:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            Chem.SanitizeMol(mol)
            Chem.MolToPDBFile(mol, str(output_path))
            print(f"Minimized molecule saved as {output_path}")
            return True
    except Exception:
        pass

    try:
        with tempfile.NamedTemporaryFile("w", suffix=".smi", delete=False) as handle:
            handle.write(smiles + "\n")
            tmp_smi = handle.name
        subprocess.run(
            [obabel, tmp_smi, "-O", str(output_path), "--gen3d", "--minimize"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        print(f"Minimized molecule saved as {output_path} via Open Babel fallback")
        return True
    except Exception:
        print(f"Energy minimization failed for SMILES: {smiles}")
        return False
    finally:
        if "tmp_smi" in locals() and os.path.exists(tmp_smi):
            os.remove(tmp_smi)


def download_pdb(pdb_id: str, download_dir: Path) -> Path:
    download_dir.mkdir(parents=True, exist_ok=True)
    pdbl = PDBList()
    pdb_file_path = pdbl.retrieve_pdb_file(pdb_id, file_format="pdb", pdir=str(download_dir))
    return Path(pdb_file_path)


def remove_hetatm(input_pdb: Path, output_pdb: Path) -> None:
    with input_pdb.open("r") as infile, output_pdb.open("w") as outfile:
        for line in infile:
            if not line.startswith("HETATM"):
                outfile.write(line)


def convert_pdb_to_pdbqt_receptor(input_pdb: Path, output_pdbqt: Path, obabel: str) -> None:
    subprocess.run(
        [obabel, "-i", "pdb", str(input_pdb), "-o", "pdbqt", "-O", str(output_pdbqt), "-xr", "-xn", "-xp"],
        check=True,
    )


def convert_pdb_to_pdbqt_ligand(input_pdb: Path, output_pdbqt: Path, obabel: str) -> None:
    subprocess.run([obabel, "-i", "pdb", str(input_pdb), "-o", "pdbqt", "-O", str(output_pdbqt), "-h"], check=True)


def run_command_with_output(command: list[str], log_file: Path) -> int:
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )
    with log_file.open("w") as log:
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            log.write(line)
            sys.stdout.flush()
    return process.wait()


def compute_pocket_box_from_residues(pdb_path: Path, residue_labels: list[str]) -> np.ndarray:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", str(pdb_path))
    residue_ids = {int(str(label).strip()) for label in residue_labels if str(label).strip().isdigit()}
    alpha_carbons = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                if residue.id[1] not in residue_ids:
                    continue
                if "CA" in residue:
                    atom = residue["CA"]
                    alpha_carbons.append([atom.coord[0], atom.coord[1], atom.coord[2]])
        break
    if not alpha_carbons:
        raise ValueError(f"No CA atoms found for P2Rank pocket residues in {pdb_path}")
    alpha_carbons = np.array(alpha_carbons)
    min_coords = np.min(alpha_carbons, axis=0)
    max_coords = np.max(alpha_carbons, axis=0)
    return max_coords - min_coords


def _default_repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def perform_docking(
    smiles_list: list[str],
    pdb_id: str,
    output_root: str | Path = "docking_results",
    config: DockingRunConfig | None = None,
) -> Path:
    """Run automated docking for a list of SMILES against a single target."""
    config = config or DockingRunConfig(output_root=Path(output_root))
    config.binary_paths.validate()

    repo_root = _default_repo_root()
    prank_executable = config.binary_paths.resolve_p2rank_executable(repo_root)

    folder_name = Path(output_root).resolve()
    receptor_name = pdb_id

    if folder_name.exists() and config.overwrite:
        shutil.rmtree(folder_name)
    folder_name.mkdir(parents=True, exist_ok=True)

    print(f"Receptor Name: {receptor_name}")
    print(f"Number of ligands: {len(smiles_list)}")

    valid_ligands: list[tuple[int, str]] = []
    for index, smiles in enumerate(smiles_list, start=1):
        pdb_filename = folder_name / f"ligand_{index}.pdb"
        if generate_minimized_pdb(smiles, pdb_filename, obabel=config.binary_paths.obabel):
            valid_ligands.append((index, smiles))
        else:
            print(f"Skipping invalid SMILES: {smiles}")

    print(f"Number of valid SMILES processed: {len(valid_ligands)}")

    downloaded_pdb_path = download_pdb(pdb_id, folder_name)
    dirty_pdb = folder_name / f"{receptor_name}_dirty.pdb"
    cleaned_pdb = folder_name / f"{receptor_name}.pdb"
    os.replace(downloaded_pdb_path, dirty_pdb)
    remove_hetatm(dirty_pdb, cleaned_pdb)
    if not config.keep_dirty_pdb:
        dirty_pdb.unlink(missing_ok=True)

    subprocess.run([str(prank_executable), "predict", "-f", str(cleaned_pdb)], check=True)

    p2rank_out = prank_executable.parent / "test_output" / f"predict_{receptor_name}"
    predictions = pd.read_csv(p2rank_out / f"{receptor_name}.pdb_predictions.csv")
    residues = pd.read_csv(p2rank_out / f"{receptor_name}.pdb_residues.csv")

    center_x = float(predictions["   center_x"].iloc[0])
    center_y = float(predictions["   center_y"].iloc[0])
    center_z = float(predictions["   center_z"].iloc[0])
    pocket1 = residues[residues[" pocket"] == 1]
    size_x, size_y, size_z = compute_pocket_box_from_residues(
        cleaned_pdb,
        pocket1[" residue_label"].tolist(),
    )

    print(center_x, center_y, center_z, size_x, size_y, size_z)

    receptor_pdbqt = folder_name / f"{receptor_name}.pdbqt"
    print("Converting receptor to PDBQT format...")
    convert_pdb_to_pdbqt_receptor(cleaned_pdb, receptor_pdbqt, config.binary_paths.obabel)
    print("Receptor conversion complete.")

    results_file = folder_name / "docking_results.txt"
    with results_file.open("w") as handle:
        handle.write("SMILES,Docking Score\n")

    for order, (index, smiles) in enumerate(valid_ligands, start=1):
        print(f"\nProcessing ligand {order} of {len(valid_ligands)}")
        print(f"SMILES: {smiles}")

        ligand_pdb = folder_name / f"ligand_{index}.pdb"
        ligand_pdbqt = folder_name / f"ligand_{index}.pdbqt"

        print("Converting ligand to PDBQT format...")
        convert_pdb_to_pdbqt_ligand(ligand_pdb, ligand_pdbqt, config.binary_paths.obabel)
        print("Ligand conversion complete.")

        output = folder_name / f"{receptor_name}_ligand_{index}.pdbqt"
        log_file = folder_name / f"vina_log_{index}.txt"
        vina_command = [
            config.binary_paths.vina,
            "--receptor",
            str(receptor_pdbqt),
            "--ligand",
            str(ligand_pdbqt),
            "--out",
            str(output),
            "--center_x",
            str(center_x),
            "--center_y",
            str(center_y),
            "--center_z",
            str(center_z),
            "--size_x",
            str(size_x),
            "--size_y",
            str(size_y),
            "--size_z",
            str(size_z),
            "--exhaustiveness",
            str(config.exhaustiveness),
            "--num_modes",
            str(config.num_modes),
        ]

        print("Starting Vina docking...")
        exit_code = run_command_with_output(vina_command, log_file)

        if exit_code == 0:
            print("Vina docking completed successfully.")
            score = "N/A"
            with log_file.open("r") as log:
                for line in log:
                    if line.startswith("   1"):
                        score = line.split()[1]
                        break
            print(f"Best docking score: {score}")
        else:
            print(f"Error running Vina for ligand {index}. Check the log file for details.")
            score = "Error"

        with results_file.open("a") as handle:
            handle.write(f"{smiles},{score}\n")

    valid_smiles_set = {smiles for _, smiles in valid_ligands}
    for smiles in smiles_list:
        if smiles not in valid_smiles_set:
            with results_file.open("a") as handle:
                handle.write(f"{smiles},Error\n")

    print(f"\nDocking complete. Results have been saved to {results_file}")
    print(f"Individual log files for each ligand are saved in the {folder_name} directory.")
    return results_file
