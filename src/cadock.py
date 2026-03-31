"""Compatibility module kept for existing scripts.

Use `dockcadd` instead of importing from `src.cadock` in new code.
"""

from dockcadd import BinaryPaths, DockingRunConfig, generate_minimized_pdb, perform_docking

__all__ = [
    "BinaryPaths",
    "DockingRunConfig",
    "generate_minimized_pdb",
    "perform_docking",
]
