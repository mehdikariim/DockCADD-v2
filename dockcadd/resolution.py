from __future__ import annotations

from typing import Iterable

import pandas as pd
import pubchempy as pcp

from .models import LigandInput


def resolve_pubchem_compounds(compounds: Iterable[dict]) -> pd.DataFrame:
    """Resolve compound names/queries to SMILES using PubChem."""
    resolved: list[dict] = []
    for item in compounds:
        name = item["compound"]
        queries = item.get("queries") or [name]
        hit = None
        used_query = None
        for query in queries:
            hits = pcp.get_compounds(query, "name")
            if hits:
                hit = hits[0]
                used_query = query
                break
        if hit is None:
            raise RuntimeError(f"Could not resolve compound via PubChem: {name}")
        smiles = getattr(hit, "smiles", None) or getattr(hit, "isomeric_smiles", None)
        if not smiles:
            raise RuntimeError(f"Resolved compound has no SMILES string: {name}")
        resolved.append(
            {
                "compound": name,
                "query_used": used_query,
                "cid": hit.cid,
                "iupac_name": getattr(hit, "iupac_name", None),
                "smiles": smiles,
            }
        )
    return pd.DataFrame(resolved)


def ligands_from_dataframe(df: pd.DataFrame) -> list[LigandInput]:
    return [LigandInput(name=row["compound"], smiles=row["smiles"]) for _, row in df.iterrows()]
