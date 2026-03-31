from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .config import BinaryPaths, DockingRunConfig
from .engine import perform_docking
from .workflows import run_matrix_from_config


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="DockCADD automated docking runner")
    subparsers = parser.add_subparsers(dest="command", required=True)

    dock_parser = subparsers.add_parser("dock", help="Dock one target against provided ligands")
    dock_parser.add_argument("--pdb-id", required=True, help="Target PDB identifier")
    dock_parser.add_argument(
        "--ligands-csv",
        type=Path,
        help="CSV file containing a SMILES column or a compound/SMILES table",
    )
    dock_parser.add_argument(
        "--smiles",
        nargs="*",
        default=[],
        help="One or more SMILES strings provided directly on the command line",
    )
    dock_parser.add_argument("--smiles-column", default="SMILES", help="Name of the SMILES column in --ligands-csv")
    dock_parser.add_argument("--output-root", type=Path, default=Path("docking_results"))
    dock_parser.add_argument("--exhaustiveness", type=int, default=8)
    dock_parser.add_argument("--num-modes", type=int, default=9)
    dock_parser.add_argument("--p2rank-dir", type=Path, help="Path to an extracted P2Rank directory")
    dock_parser.add_argument("--vina", default="vina", help="Vina executable name or path")
    dock_parser.add_argument("--obabel", default="obabel", help="Open Babel executable name or path")

    matrix_parser = subparsers.add_parser("matrix", help="Run a config-driven ligand/target matrix")
    matrix_parser.add_argument("--config-json", type=Path, required=True, help="JSON config file for a matrix run")
    return parser


def _load_smiles(args: argparse.Namespace) -> list[str]:
    smiles = list(args.smiles)
    if args.ligands_csv:
        df = pd.read_csv(args.ligands_csv)
        if args.smiles_column not in df.columns:
            raise ValueError(f"Column '{args.smiles_column}' was not found in {args.ligands_csv}")
        smiles.extend(df[args.smiles_column].dropna().astype(str).tolist())
    if not smiles:
        raise ValueError("Provide ligands with --smiles and/or --ligands-csv.")
    return smiles


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "dock":
        smiles = _load_smiles(args)
        config = DockingRunConfig(
            output_root=args.output_root,
            exhaustiveness=args.exhaustiveness,
            num_modes=args.num_modes,
            binary_paths=BinaryPaths(
                vina=args.vina,
                obabel=args.obabel,
                p2rank_dir=args.p2rank_dir,
            ),
        )
        perform_docking(smiles, args.pdb_id, output_root=args.output_root, config=config)
    elif args.command == "matrix":
        run_matrix_from_config(args.config_json)
    else:
        raise ValueError(f"Unsupported command: {args.command}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
