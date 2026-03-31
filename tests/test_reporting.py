from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


class ReportingTests(unittest.TestCase):
    def test_write_report_bundle_creates_summary_sheet(self) -> None:
        from dockcadd.reporting import write_report_bundle

        df = pd.DataFrame(
            [
                {"Target": "A", "PDB ID": "1ABC", "Compound": "Lig1", "SMILES": "CCO", "Binding affinity (kcal/mol)": -7.1},
                {"Target": "A", "PDB ID": "1ABC", "Compound": "Lig2", "SMILES": "CCN", "Binding affinity (kcal/mol)": -6.2},
                {"Target": "B", "PDB ID": "2XYZ", "Compound": "Lig3", "SMILES": "CCC", "Binding affinity (kcal/mol)": -8.0},
            ]
        )
        resolved = pd.DataFrame([{"compound": "Lig1", "smiles": "CCO"}])
        with tempfile.TemporaryDirectory() as tmpdir:
            xlsx = Path(tmpdir) / "report.xlsx"
            csv = Path(tmpdir) / "report.csv"
            write_report_bundle(df, xlsx, csv, resolved_df=resolved)
            self.assertTrue(xlsx.exists())
            self.assertTrue(csv.exists())
            with pd.ExcelFile(xlsx) as xl:
                self.assertIn("DockingResults", xl.sheet_names)
                self.assertIn("TargetSummary", xl.sheet_names)
                self.assertIn("CompoundResolution", xl.sheet_names)


if __name__ == "__main__":
    unittest.main()
