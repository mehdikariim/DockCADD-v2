from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .config import BinaryPaths, DockingRunConfig
from .engine import (
    load_smiles_from_pdb_files,
    load_smiles_from_sdf_files,
    perform_docking,
    perform_redocking,
)
from .workflows import run_matrix_from_config


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="DockCADD automated docking runner")
    subparsers = parser.add_subparsers(dest="command", required=True)

    dock_parser = subparsers.add_parser("dock", help="Dock one target against provided ligands")
    _add_target_arguments(dock_parser)
    _add_common_run_arguments(dock_parser)
    dock_parser.add_argument(
        "--ligands-csv",
        type=Path,
        help="CSV file containing a SMILES column.",
    )
    dock_parser.add_argument(
        "--ligands-sdf",
        type=Path,
        nargs="*",
        default=[],
        help="One or more SDF files. All molecules found in each file are docked.",
    )
    dock_parser.add_argument(
        "--ligands-pdb",
        type=Path,
        nargs="*",
        default=[],
        help="One or more ligand PDB files. SMILES are derived automatically.",
    )
    dock_parser.add_argument(
        "--smiles",
        nargs="*",
        default=[],
        help="One or more SMILES strings provided directly on the command line.",
    )
    dock_parser.add_argument(
        "--smiles-column",
        default="SMILES",
        help="Name of the SMILES column in --ligands-csv.",
    )

    redock_parser = subparsers.add_parser("redock", help="Redock a known ligand into a target")
    _add_target_arguments(redock_parser)
    _add_common_run_arguments(redock_parser)
    redock_parser.add_argument(
        "--ligand-file",
        type=Path,
        required=True,
        help="Ligand structure file to redock (.pdb or .sdf).",
    )

    matrix_parser = subparsers.add_parser("matrix", help="Run a config-driven ligand/target matrix")
    matrix_parser.add_argument("--config-json", type=Path, required=True, help="JSON config file for a matrix run")
    return parser


def _add_target_arguments(parser: argparse.ArgumentParser) -> None:
    target_group = parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument("--pdb-id", help="Target PDB identifier")
    target_group.add_argument("--receptor-pdb", type=Path, help="Local receptor PDB file")


def _add_common_run_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--output-root", type=Path, default=Path("docking_results"))
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--num-modes", type=int, default=9)
    parser.add_argument("--disable-pdbfixer", action="store_true", help="Skip PDBFixer cleanup if installed")
    parser.add_argument("--keep-dirty-pdb", action="store_true", help="Keep the unprocessed receptor PDB copy")
    parser.add_argument("--no-overwrite", action="store_true", help="Fail instead of replacing an existing output folder")
    parser.add_argument("--p2rank-dir", type=Path, help="Path to an extracted P2Rank directory")
    parser.add_argument("--vina", default="vina", help="Vina executable name or path")
    parser.add_argument("--obabel", default="obabel", help="Open Babel executable name or path")


def _load_smiles(args: argparse.Namespace) -> list[str]:
    smiles = list(args.smiles)

    if args.ligands_csv:
        df = pd.read_csv(args.ligands_csv)
        if args.smiles_column not in df.columns:
            raise ValueError(f"Column '{args.smiles_column}' was not found in {args.ligands_csv}")
        smiles.extend(df[args.smiles_column].dropna().astype(str).tolist())

    if args.ligands_sdf:
        smiles.extend(load_smiles_from_sdf_files(args.ligands_sdf))

    if args.ligands_pdb:
        smiles.extend(load_smiles_from_pdb_files(args.ligands_pdb, args.obabel))

    if not smiles:
        raise ValueError("Provide ligands with --smiles, --ligands-csv, --ligands-sdf, and/or --ligands-pdb.")
    return smiles


def _build_config(args: argparse.Namespace) -> DockingRunConfig:
    return DockingRunConfig(
        output_root=args.output_root,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        use_pdbfixer=not args.disable_pdbfixer,
        keep_dirty_pdb=args.keep_dirty_pdb,
        overwrite=not args.no_overwrite,
        binary_paths=BinaryPaths(
            vina=args.vina,
            obabel=args.obabel,
            p2rank_dir=args.p2rank_dir,
        ),
    )


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "dock":
        smiles = _load_smiles(args)
        config = _build_config(args)
        perform_docking(
            smiles,
            args.pdb_id,
            receptor_pdb=args.receptor_pdb,
            output_root=args.output_root,
            config=config,
        )
    elif args.command == "redock":
        config = _build_config(args)
        perform_redocking(
            args.ligand_file,
            args.pdb_id,
            receptor_pdb=args.receptor_pdb,
            output_root=args.output_root,
            config=config,
        )
    elif args.command == "matrix":
        run_matrix_from_config(args.config_json)
    else:
        raise ValueError(f"Unsupported command: {args.command}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
