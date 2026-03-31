from __future__ import annotations

from pathlib import Path
import json

import pandas as pd

from .config import BinaryPaths, DockingRunConfig
from .engine import perform_docking
from .reporting import write_report_bundle
from .resolution import resolve_pubchem_compounds


def _target_run_id(target: dict) -> str:
    if target.get("pdb_id"):
        return str(target["pdb_id"])
    if target.get("receptor_pdb"):
        return Path(target["receptor_pdb"]).stem
    raise ValueError("Each target must define either 'pdb_id' or 'receptor_pdb'.")


def assemble_matrix_results(
    resolved_df: pd.DataFrame,
    targets: list[dict],
    output_root: Path,
) -> pd.DataFrame:
    name_by_smiles = dict(zip(resolved_df["smiles"], resolved_df["compound"]))
    rows: list[pd.DataFrame] = []
    for target in targets:
        target_id = _target_run_id(target)
        result_file = output_root / target_id / "docking_results.txt"
        df = pd.read_csv(result_file)
        df["Compound"] = df["SMILES"].map(name_by_smiles)
        df["Target"] = target["target"]
        df["PDB ID"] = target.get("pdb_id", target_id)
        df.rename(columns={"Docking Score": "Binding affinity (kcal/mol)"}, inplace=True)
        rows.append(df[["Target", "PDB ID", "Compound", "SMILES", "Binding affinity (kcal/mol)"]])
    final_df = pd.concat(rows, ignore_index=True)
    final_df["Binding affinity (kcal/mol)"] = pd.to_numeric(
        final_df["Binding affinity (kcal/mol)"], errors="coerce"
    )
    final_df.sort_values(by=["Target", "Binding affinity (kcal/mol)"], inplace=True)
    return final_df


def run_matrix_from_resolution(
    resolved_df: pd.DataFrame,
    targets: list[dict],
    output_root: Path,
    config: DockingRunConfig | None = None,
) -> pd.DataFrame:
    output_root.mkdir(parents=True, exist_ok=True)
    smiles_list = resolved_df["smiles"].tolist()

    for target in targets:
        target_dir = output_root / _target_run_id(target)
        run_config = config or DockingRunConfig(output_root=target_dir)
        perform_docking(
            smiles_list,
            target.get("pdb_id"),
            receptor_pdb=target.get("receptor_pdb"),
            output_root=target_dir,
            config=run_config,
        )
    return assemble_matrix_results(resolved_df=resolved_df, targets=targets, output_root=output_root)


def run_matrix_from_config(config_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    config_data = json.loads(config_path.read_text(encoding="utf-8"))
    base_dir = config_path.resolve().parent
    output_root = _resolve_config_path(base_dir, config_data["output_root"])
    report_xlsx = _resolve_config_path(base_dir, config_data["report_xlsx"])
    report_csv = _resolve_config_path(base_dir, config_data["report_csv"])
    dock_config = DockingRunConfig(
        output_root=output_root,
        exhaustiveness=config_data.get("exhaustiveness", 8),
        num_modes=config_data.get("num_modes", 9),
        use_pdbfixer=config_data.get("use_pdbfixer", True),
        keep_dirty_pdb=config_data.get("keep_dirty_pdb", True),
        overwrite=config_data.get("overwrite", True),
        binary_paths=BinaryPaths(
            vina=config_data.get("vina", "vina"),
            obabel=config_data.get("obabel", "obabel"),
            p2rank_dir=_resolve_config_path(base_dir, config_data["p2rank_dir"]) if config_data.get("p2rank_dir") else None,
        ),
    )
    resolved_df = resolve_pubchem_compounds(config_data["compounds"])
    matrix_df = run_matrix_from_resolution(
        resolved_df=resolved_df,
        targets=config_data["targets"],
        output_root=output_root,
        config=dock_config,
    )
    write_report_bundle(matrix_df, report_xlsx, report_csv, resolved_df=resolved_df)
    return matrix_df, resolved_df


def _resolve_config_path(base_dir: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (base_dir / path).resolve()
