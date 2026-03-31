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

try:
    from pdbfixer import PDBFixer
    try:
        from openmm.app import PDBFile
    except ImportError:  # pragma: no cover
        from simtk.openmm.app import PDBFile  # type: ignore
except ImportError:  # pragma: no cover
    PDBFixer = None
    PDBFile = None


def generate_minimized_pdb(smiles: str, pdb_filename: str | Path, *, obabel: str = "obabel") -> bool:
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


def fix_receptor_structure(input_pdb: Path, output_pdb: Path, *, use_pdbfixer: bool = True) -> None:
    if use_pdbfixer and PDBFixer is not None and PDBFile is not None:
        try:
            fixer = PDBFixer(filename=str(input_pdb))
            fixer.removeHeterogens(keepWater=False)
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(7.0)
            with output_pdb.open("w") as handle:
                PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
            return
        except Exception as exc:
            print(f"PDBFixer cleanup failed, falling back to HETATM filtering: {exc}")
    remove_hetatm(input_pdb, output_pdb)


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


def _prepare_receptor(
    *,
    folder_name: Path,
    config: DockingRunConfig,
    pdb_id: str | None = None,
    receptor_pdb: str | Path | None = None,
) -> tuple[str, Path]:
    if bool(pdb_id) == bool(receptor_pdb):
        raise ValueError("Provide exactly one of pdb_id or receptor_pdb.")

    if receptor_pdb is not None:
        receptor_source = Path(receptor_pdb).resolve()
        receptor_name = receptor_source.stem
        dirty_pdb = folder_name / f"{receptor_name}_dirty.pdb"
        shutil.copy2(receptor_source, dirty_pdb)
    else:
        receptor_name = str(pdb_id)
        downloaded_pdb_path = download_pdb(receptor_name, folder_name)
        dirty_pdb = folder_name / f"{receptor_name}_dirty.pdb"
        os.replace(downloaded_pdb_path, dirty_pdb)

    cleaned_pdb = folder_name / f"{receptor_name}.pdb"
    fix_receptor_structure(dirty_pdb, cleaned_pdb, use_pdbfixer=config.use_pdbfixer)
    if not config.keep_dirty_pdb:
        dirty_pdb.unlink(missing_ok=True)
    return receptor_name, cleaned_pdb


def load_smiles_from_sdf_files(sdf_files: list[Path]) -> list[str]:
    smiles: list[str] = []
    for sdf_file in sdf_files:
        supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False)
        for mol in supplier:
            if mol is None:
                continue
            smiles.append(Chem.MolToSmiles(Chem.RemoveHs(mol)))
    return smiles


def smiles_from_structure_file(path: Path, obabel: str) -> str:
    suffix = path.suffix.lower()
    if suffix == ".pdb":
        mol = Chem.MolFromPDBFile(str(path), removeHs=False)
        if mol is not None:
            return Chem.MolToSmiles(Chem.RemoveHs(mol))
    elif suffix == ".sdf":
        supplier = Chem.SDMolSupplier(str(path), removeHs=False)
        mol = next((item for item in supplier if item is not None), None)
        if mol is not None:
            return Chem.MolToSmiles(Chem.RemoveHs(mol))

    result = subprocess.run(
        [obabel, str(path), "-ocan"],
        check=True,
        capture_output=True,
        text=True,
    )
    line = next((line.strip() for line in result.stdout.splitlines() if line.strip()), "")
    if not line:
        raise RuntimeError(f"Could not derive SMILES from {path}")
    return line.split()[0]


def load_smiles_from_pdb_files(pdb_files: list[Path], obabel: str) -> list[str]:
    return [smiles_from_structure_file(path, obabel) for path in pdb_files]


def convert_structure_file_to_pdb(input_path: Path, output_pdb: Path, obabel: str) -> None:
    suffix = input_path.suffix.lower()
    if suffix == ".pdb":
        shutil.copy2(input_path, output_pdb)
        return
    if suffix == ".sdf":
        subprocess.run([obabel, str(input_path), "-O", str(output_pdb)], check=True)
        return
    raise ValueError(f"Unsupported ligand file for redocking: {input_path}")


def _predict_box(cleaned_pdb: Path, receptor_name: str, prank_executable: Path) -> tuple[float, float, float, float, float, float]:
    subprocess.run([str(prank_executable), "predict", "-f", str(cleaned_pdb)], check=True)
    p2rank_out = prank_executable.parent / "test_output" / f"predict_{receptor_name}"
    predictions = pd.read_csv(p2rank_out / f"{receptor_name}.pdb_predictions.csv")
    residues = pd.read_csv(p2rank_out / f"{receptor_name}.pdb_residues.csv")
    center_x = float(predictions["   center_x"].iloc[0])
    center_y = float(predictions["   center_y"].iloc[0])
    center_z = float(predictions["   center_z"].iloc[0])
    pocket1 = residues[residues[" pocket"] == 1]
    size_x, size_y, size_z = compute_pocket_box_from_residues(cleaned_pdb, pocket1[" residue_label"].tolist())
    return center_x, center_y, center_z, float(size_x), float(size_y), float(size_z)


def _dock_prepared_ligands(
    *,
    folder_name: Path,
    receptor_name: str,
    receptor_pdbqt: Path,
    ligand_entries: list[tuple[int, str]],
    center_x: float,
    center_y: float,
    center_z: float,
    size_x: float,
    size_y: float,
    size_z: float,
    config: DockingRunConfig,
) -> Path:
    results_file = folder_name / "docking_results.txt"
    with results_file.open("w") as handle:
        handle.write("SMILES,Docking Score\n")

    for order, (index, smiles) in enumerate(ligand_entries, start=1):
        print(f"\nProcessing ligand {order} of {len(ligand_entries)}")
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
            "--receptor", str(receptor_pdbqt),
            "--ligand", str(ligand_pdbqt),
            "--out", str(output),
            "--center_x", str(center_x),
            "--center_y", str(center_y),
            "--center_z", str(center_z),
            "--size_x", str(size_x),
            "--size_y", str(size_y),
            "--size_z", str(size_z),
            "--exhaustiveness", str(config.exhaustiveness),
            "--num_modes", str(config.num_modes),
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
        else:
            score = "Error"
        with results_file.open("a") as handle:
            handle.write(f"{smiles},{score}\n")
    return results_file


def perform_docking(
    smiles_list: list[str],
    pdb_id: str | None = None,
    *,
    receptor_pdb: str | Path | None = None,
    output_root: str | Path = "docking_results",
    config: DockingRunConfig | None = None,
) -> Path:
    config = config or DockingRunConfig(output_root=Path(output_root))
    config.binary_paths.validate()

    repo_root = _default_repo_root()
    prank_executable = config.binary_paths.resolve_p2rank_executable(repo_root)
    folder_name = Path(output_root).resolve()

    if folder_name.exists() and config.overwrite:
        shutil.rmtree(folder_name)
    folder_name.mkdir(parents=True, exist_ok=True)

    receptor_name, cleaned_pdb = _prepare_receptor(
        folder_name=folder_name,
        config=config,
        pdb_id=pdb_id,
        receptor_pdb=receptor_pdb,
    )
    print(f"Receptor Name: {receptor_name}")
    print(f"Number of ligands: {len(smiles_list)}")

    valid_ligands: list[tuple[int, str]] = []
    for index, smiles in enumerate(smiles_list, start=1):
        pdb_filename = folder_name / f"ligand_{index}.pdb"
        if generate_minimized_pdb(smiles, pdb_filename, obabel=config.binary_paths.obabel):
            valid_ligands.append((index, smiles))
        else:
            print(f"Skipping invalid SMILES: {smiles}")

    center_x, center_y, center_z, size_x, size_y, size_z = _predict_box(cleaned_pdb, receptor_name, prank_executable)
    print(center_x, center_y, center_z, size_x, size_y, size_z)

    receptor_pdbqt = folder_name / f"{receptor_name}.pdbqt"
    print("Converting receptor to PDBQT format...")
    convert_pdb_to_pdbqt_receptor(cleaned_pdb, receptor_pdbqt, config.binary_paths.obabel)
    print("Receptor conversion complete.")

    results_file = _dock_prepared_ligands(
        folder_name=folder_name,
        receptor_name=receptor_name,
        receptor_pdbqt=receptor_pdbqt,
        ligand_entries=valid_ligands,
        center_x=center_x,
        center_y=center_y,
        center_z=center_z,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        config=config,
    )

    valid_smiles_set = {smiles for _, smiles in valid_ligands}
    for smiles in smiles_list:
        if smiles not in valid_smiles_set:
            with results_file.open("a") as handle:
                handle.write(f"{smiles},Error\n")

    print(f"\nDocking complete. Results have been saved to {results_file}")
    return results_file


def perform_redocking(
    ligand_file: str | Path,
    pdb_id: str | None = None,
    *,
    receptor_pdb: str | Path | None = None,
    output_root: str | Path = "redocking_results",
    config: DockingRunConfig | None = None,
) -> Path:
    config = config or DockingRunConfig(output_root=Path(output_root))
    config.binary_paths.validate()

    repo_root = _default_repo_root()
    prank_executable = config.binary_paths.resolve_p2rank_executable(repo_root)
    folder_name = Path(output_root).resolve()
    if folder_name.exists() and config.overwrite:
        shutil.rmtree(folder_name)
    folder_name.mkdir(parents=True, exist_ok=True)

    receptor_name, cleaned_pdb = _prepare_receptor(
        folder_name=folder_name,
        config=config,
        pdb_id=pdb_id,
        receptor_pdb=receptor_pdb,
    )
    ligand_path = Path(ligand_file).resolve()
    ligand_pdb = folder_name / "ligand_1.pdb"
    convert_structure_file_to_pdb(ligand_path, ligand_pdb, config.binary_paths.obabel)
    ligand_smiles = smiles_from_structure_file(ligand_path, config.binary_paths.obabel)

    center_x, center_y, center_z, size_x, size_y, size_z = _predict_box(cleaned_pdb, receptor_name, prank_executable)
    receptor_pdbqt = folder_name / f"{receptor_name}.pdbqt"
    convert_pdb_to_pdbqt_receptor(cleaned_pdb, receptor_pdbqt, config.binary_paths.obabel)

    results_file = _dock_prepared_ligands(
        folder_name=folder_name,
        receptor_name=receptor_name,
        receptor_pdbqt=receptor_pdbqt,
        ligand_entries=[(1, ligand_smiles)],
        center_x=center_x,
        center_y=center_y,
        center_z=center_z,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        config=config,
    )
    print(f"\nRedocking complete. Results have been saved to {results_file}")
    return results_file
