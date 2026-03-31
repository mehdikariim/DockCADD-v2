from pathlib import Path

import pandas as pd

from dockcadd import perform_docking


COMPOUND = "4-hydroxy-beta-thujone"
SMILES = "CC(C)[C@@]12C[C@H]1[C@@](C(=O)C2)(C)O"
TARGETS = [
    {"target": "Antifungal", "pdb_id": "5TZ1"},
    {"target": "Antioxidant", "pdb_id": "1US0"},
    {"target": "Larvicidal", "pdb_id": "5W1U"},
]


def main() -> None:
    desktop_dir = Path("/mnt/c/Users/uid-1180/Desktop")
    output_root = desktop_dir / "chebbak_dockcaddv2_results_patch"
    output_root.mkdir(parents=True, exist_ok=True)

    patched = []
    for target in TARGETS:
        target_dir = output_root / target["pdb_id"]
        perform_docking([SMILES], target["pdb_id"], output_root=str(target_dir))
        df = pd.read_csv(target_dir / "docking_results.txt")
        score = df.iloc[0]["Docking Score"]
        patched.append(
            {
                "Target": target["target"],
                "PDB ID": target["pdb_id"],
                "Compound": COMPOUND,
                "SMILES": SMILES,
                "Binding affinity (kcal/mol)": score,
            }
        )

    patched_df = pd.DataFrame(patched)
    csv_path = desktop_dir / "chebbak_dockcaddv2_all_compounds.csv"
    xlsx_path = desktop_dir / "chebbak_dockcaddv2_all_compounds.xlsx"

    main_df = pd.read_csv(csv_path, dtype=str)
    for _, row in patched_df.iterrows():
        mask = (main_df["PDB ID"] == row["PDB ID"]) & (main_df["Compound"] == row["Compound"])
        main_df.loc[mask, "Binding affinity (kcal/mol)"] = f"{float(row['Binding affinity (kcal/mol)']):.3f}"
    main_df.to_csv(csv_path, index=False)

    resolution_df = pd.read_excel(xlsx_path, sheet_name="CompoundResolution")
    docking_df = main_df.copy()
    docking_df["Binding affinity (kcal/mol)"] = pd.to_numeric(
        docking_df["Binding affinity (kcal/mol)"], errors="coerce"
    )
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        docking_df.to_excel(writer, index=False, sheet_name="DockingResults")
        resolution_df.to_excel(writer, index=False, sheet_name="CompoundResolution")

    print(csv_path)
    print(xlsx_path)


if __name__ == "__main__":
    main()
