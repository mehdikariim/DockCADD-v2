from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import shutil
import sys


@dataclass(slots=True)
class BinaryPaths:
    """External tool locations used by DockCADD."""

    vina: str = "vina"
    obabel: str = "obabel"
    p2rank_dir: Path | None = None

    def resolve_p2rank_executable(self, repo_root: Path) -> Path:
        base = self.p2rank_dir or repo_root / "p2rank_2.4.2"
        executable = "prank.bat" if sys.platform.startswith("win") else "prank"
        prank = base / executable
        if not prank.exists():
            raise FileNotFoundError(
                f"P2Rank executable was not found at {prank}. "
                "Install P2Rank or point BinaryPaths.p2rank_dir to the extracted directory."
            )
        return prank

    def validate(self) -> None:
        missing = [tool for tool in (self.vina, self.obabel) if shutil.which(tool) is None]
        if missing:
            raise FileNotFoundError(
                "Missing required executables in PATH: " + ", ".join(missing)
            )


@dataclass(slots=True)
class DockingRunConfig:
    """User-facing knobs for a docking run."""

    output_root: Path = field(default_factory=lambda: Path("docking_results"))
    exhaustiveness: int = 8
    num_modes: int = 9
    use_pdbfixer: bool = True
    keep_dirty_pdb: bool = True
    overwrite: bool = True
    binary_paths: BinaryPaths = field(default_factory=BinaryPaths)
