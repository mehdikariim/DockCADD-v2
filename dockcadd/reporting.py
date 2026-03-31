from __future__ import annotations

from pathlib import Path
from datetime import datetime

import pandas as pd


def normalize_results_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    data = df.copy()
    data["Binding affinity (kcal/mol)"] = pd.to_numeric(
        data["Binding affinity (kcal/mol)"], errors="coerce"
    )
    return data


def build_target_summary(df: pd.DataFrame, top_n: int = 3) -> pd.DataFrame:
    summary_rows: list[dict] = []
    for (target, pdb_id), group in df.groupby(["Target", "PDB ID"], sort=True):
        ordered = group.sort_values("Binding affinity (kcal/mol)").head(top_n)
        row = {
            "Target": target,
            "PDB ID": pdb_id,
            "Top hits": "; ".join(
                f"{compound} ({score:.3f})"
                for compound, score in zip(
                    ordered["Compound"],
                    ordered["Binding affinity (kcal/mol)"],
                )
                if pd.notna(score)
            ),
        }
        summary_rows.append(row)
    return pd.DataFrame(summary_rows)


def write_report_bundle(
    docking_df: pd.DataFrame,
    output_xlsx: Path,
    output_csv: Path,
    *,
    resolved_df: pd.DataFrame | None = None,
    top_n: int = 3,
) -> tuple[Path, Path]:
    docking_df = normalize_results_dataframe(docking_df)
    summary_df = build_target_summary(docking_df, top_n=top_n)
    output_xlsx.parent.mkdir(parents=True, exist_ok=True)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    csv_target = _safe_output_target(output_csv)
    xlsx_target = _safe_output_target(output_xlsx)
    docking_df.to_csv(csv_target, index=False)
    with pd.ExcelWriter(xlsx_target, engine="openpyxl") as writer:
        docking_df.to_excel(writer, index=False, sheet_name="DockingResults")
        summary_df.to_excel(writer, index=False, sheet_name="TargetSummary")
        if resolved_df is not None:
            resolved_df.to_excel(writer, index=False, sheet_name="CompoundResolution")
    return xlsx_target, csv_target


def _safe_output_target(path: Path) -> Path:
    try:
        with path.open("ab"):
            return path
    except PermissionError:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return path.with_name(f"{path.stem}_{stamp}{path.suffix}")
