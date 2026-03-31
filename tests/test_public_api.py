from __future__ import annotations

import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


class PublicApiTests(unittest.TestCase):
    def test_package_exports(self) -> None:
        import dockcadd

        self.assertIn("perform_docking", dockcadd.__all__)
        self.assertIn("perform_redocking", dockcadd.__all__)
        self.assertIn("load_smiles_from_sdf_files", dockcadd.__all__)
        self.assertTrue(callable(dockcadd.perform_docking))
        self.assertTrue(callable(dockcadd.perform_redocking))

    def test_dock_cli_parser_supports_local_receptor_and_sdf(self) -> None:
        from dockcadd.cli import build_parser

        parser = build_parser()
        args = parser.parse_args(
            [
                "dock",
                "--receptor-pdb",
                "receptor.pdb",
                "--ligands-sdf",
                "ligands.sdf",
                "--output-root",
                "outputs/test",
            ]
        )
        self.assertEqual(args.command, "dock")
        self.assertEqual(args.receptor_pdb, Path("receptor.pdb"))
        self.assertEqual(args.ligands_sdf, [Path("ligands.sdf")])

    def test_redock_cli_parser_exists(self) -> None:
        from dockcadd.cli import build_parser

        parser = build_parser()
        args = parser.parse_args(
            [
                "redock",
                "--pdb-id",
                "5TZ1",
                "--ligand-file",
                "reference.sdf",
            ]
        )
        self.assertEqual(args.command, "redock")
        self.assertEqual(args.pdb_id, "5TZ1")
        self.assertEqual(args.ligand_file, Path("reference.sdf"))

    def test_package_exports_matrix_helpers(self) -> None:
        import dockcadd

        self.assertIn("run_matrix_from_config", dockcadd.__all__)
        self.assertIn("write_report_bundle", dockcadd.__all__)


if __name__ == "__main__":
    unittest.main()
