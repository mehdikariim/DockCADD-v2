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
        self.assertTrue(callable(dockcadd.perform_docking))

    def test_legacy_wrappers_forward_to_new_package(self) -> None:
        from dockcadd import perform_docking as new_impl
        from src.cadock import perform_docking as src_impl
        from caddock.docking import perform_docking as legacy_impl

        self.assertIs(new_impl, src_impl)
        self.assertIs(new_impl, legacy_impl)

    def test_cli_parser_exists(self) -> None:
        from dockcadd.cli import build_parser

        parser = build_parser()
        args = parser.parse_args(["dock", "--pdb-id", "5TZ1", "--smiles", "CCO"])
        self.assertEqual(args.command, "dock")
        self.assertEqual(args.pdb_id, "5TZ1")
        self.assertEqual(args.smiles, ["CCO"])

    def test_package_exports_matrix_helpers(self) -> None:
        import dockcadd

        self.assertIn("run_matrix_from_config", dockcadd.__all__)
        self.assertIn("write_report_bundle", dockcadd.__all__)


if __name__ == "__main__":
    unittest.main()
