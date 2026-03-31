from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class LigandInput:
    name: str
    smiles: str


@dataclass(slots=True)
class TargetInput:
    name: str
    pdb_id: str


@dataclass(slots=True)
class MatrixRow:
    target: str
    pdb_id: str
    compound: str
    smiles: str
    binding_affinity_kcal_mol: str
