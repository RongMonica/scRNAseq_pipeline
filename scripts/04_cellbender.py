"""Run CellBender on raw 10x-style H5 matrices.

Input:
    data/raw_h5/<sample>/raw_feature_bc_matrix.h5

Output:
    data/clean_data/<sample>/<sample>_cellbender.h5

Usage:
    python scripts/04_cellbender.py

Pass CellBender options with CELLBENDER_EXTRA_ARGS:
    CELLBENDER_EXTRA_ARGS="--expected-cells 5000 --total-droplets-included 20000" \
        python scripts/04_cellbender.py
"""

from pathlib import Path
import os
import shlex
import shutil
import subprocess

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
CLEAN_DATA_DIR = DATA_DIR / "clean_data"
RAW_H5_DIR = DATA_DIR / "raw_h5"
PREPROCESS_DIR = REPO_ROOT / "results" / "preprocess"


def sample_id_from_h5_path(path: Path) -> str:
    # STARsolo paths are usually:
    # results/preprocess/<sample>/Solo.out/Gene/raw/raw_feature_bc_matrix.h5
    if (
        path.parent.name == "raw"
        and len(path.parents) > 4
        and path.parents[4] == PREPROCESS_DIR
    ):
        return path.parents[3].name

    return path.parent.name


def find_raw_h5_files() -> dict[str, Path]:
    sample_files: dict[str, Path] = {}
    search_roots = [RAW_H5_DIR, PREPROCESS_DIR, DATA_DIR]

    for root in search_roots:
        if not root.exists():
            continue
        for path in sorted(root.rglob("raw_feature_bc_matrix.h5")):
            sample_files.setdefault(sample_id_from_h5_path(path), path)
    return sample_files


def run_cellbender(sample_id: str, raw_h5_path: Path) -> Path:
    output_path = CLEAN_DATA_DIR / sample_id / f"{sample_id}_cellbender.h5"

    if output_path.exists():
        print(f"Skipping {sample_id}; found existing {output_path}")
        return output_path

    output_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "cellbender",
        "remove-background",
        "--input",
        str(raw_h5_path),
        "--output",
        str(output_path),
    ]
    cmd.extend(shlex.split(os.environ.get("CELLBENDER_EXTRA_ARGS", "")))

    print(f"Running CellBender for {sample_id}")
    subprocess.run(cmd, check=True)
    return output_path


def main() -> None:
    raw_files = find_raw_h5_files()
    if not raw_files:
        raise FileNotFoundError("No raw_feature_bc_matrix.h5 files found")

    for sample_id, raw_h5_path in raw_files.items():
        run_cellbender(sample_id, raw_h5_path)


if __name__ == "__main__":
    main()
