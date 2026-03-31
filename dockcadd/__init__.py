"""Public package entrypoints for DockCADD."""

from .config import BinaryPaths, DockingRunConfig
from .engine import generate_minimized_pdb, perform_docking
from .reporting import write_report_bundle
from .resolution import resolve_pubchem_compounds
from .workflows import assemble_matrix_results, run_matrix_from_config, run_matrix_from_resolution

__all__ = [
    "BinaryPaths",
    "DockingRunConfig",
    "generate_minimized_pdb",
    "perform_docking",
    "resolve_pubchem_compounds",
    "assemble_matrix_results",
    "run_matrix_from_config",
    "run_matrix_from_resolution",
    "write_report_bundle",
]
