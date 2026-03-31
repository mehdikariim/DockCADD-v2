from __future__ import annotations

from pathlib import Path

from dockcadd import run_matrix_from_config


def main() -> None:
    config_path = Path(__file__).resolve().parent / "examples" / "chebbak_matrix.json"
    matrix_df, resolved_df = run_matrix_from_config(config_path)
    print(config_path)
    print(len(matrix_df))
    print(len(resolved_df))


if __name__ == "__main__":
    main()
