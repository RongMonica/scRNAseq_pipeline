"""Convert arranged 10x matrix/barcode/feature triplet files into H5 files.

Input:
    data/raw/<sample>/matrix.mtx[.gz]
    data/raw/<sample>/barcodes.tsv[.gz]
    data/raw/<sample>/features.tsv[.gz]
    data/processed/<sample>/matrix.mtx[.gz]
    data/processed/<sample>/barcodes.tsv[.gz]
    data/processed/<sample>/features.tsv[.gz]

Output:
    data/raw_h5/<sample>/raw_feature_bc_matrix.h5
    data/filtered_h5/<sample>/filtered_feature_bc_matrix.h5

Usage:
    python3 scripts/02_make_raw_h5.py
    SAMPLE_NAME=<sample> python3 scripts/02_make_raw_h5.py
    COUNT_TYPE=raw python3 scripts/02_make_raw_h5.py
    COUNT_TYPE=processed python3 scripts/02_make_raw_h5.py
"""

import gzip
import os
from pathlib import Path

import h5py
import numpy as np
from scipy.io import mmread
from scipy.sparse import csc_matrix


REPO_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = REPO_ROOT / "data" / "raw"
RAW_H5_DIR = REPO_ROOT / "data" / "raw_h5"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
FILTERED_H5_DIR = REPO_ROOT / "data" / "filtered_h5"

COUNT_CONFIGS = {
    "raw": {
        "input_dir": RAW_DIR,
        "output_dir": RAW_H5_DIR,
        "h5_name": "raw_feature_bc_matrix.h5",
    },
    "processed": {
        "input_dir": PROCESSED_DIR,
        "output_dir": FILTERED_H5_DIR,
        "h5_name": "filtered_feature_bc_matrix.h5",
    },
}


def first_existing(*paths: Path) -> Path:
    """Return the first path that exists, allowing .gz or uncompressed inputs."""
    for path in paths:
        if path.exists():
            return path
    raise FileNotFoundError(f"None of these files exist: {', '.join(map(str, paths))}")


def read_lines(path: Path) -> list[str]:
    """Read a text file that may be plain text or gzip-compressed."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle]


def is_10x_h5_valid(path: Path) -> bool:
    """Check whether a path contains a readable 10x-style H5 file."""
    if not path.is_file():
        return False

    try:
        with h5py.File(path, "r") as h5:
            if "matrix" not in h5:
                return False
            matrix_group = h5["matrix"]

            required_matrix_datasets = {"barcodes", "data", "indices", "indptr", "shape"}
            if not required_matrix_datasets.issubset(matrix_group.keys()):
                return False

            if "features" not in matrix_group:
                return False
            features_group = matrix_group["features"]

            required_feature_datasets = {"id", "name", "feature_type", "genome", "_all_tag_keys"}
            if not required_feature_datasets.issubset(features_group.keys()):
                return False

            shape = matrix_group["shape"][()]
            if len(shape) != 2:
                return False

            n_features, n_barcodes = shape
            if len(matrix_group["barcodes"]) != n_barcodes:
                return False
            if len(features_group["id"]) != n_features:
                return False
            if len(features_group["name"]) != n_features:
                return False
            if len(features_group["feature_type"]) != n_features:
                return False
            if len(matrix_group["data"]) != len(matrix_group["indices"]):
                return False
            if len(matrix_group["indptr"]) != n_barcodes + 1:
                return False
        return True
    except (OSError, KeyError, TypeError, ValueError):
        return False


def write_h5(sample_dir: Path, output_path: Path) -> None:
    # Find the three required 10x-style files for this sample.
    matrix_path = first_existing(sample_dir / "matrix.mtx.gz", sample_dir / "matrix.mtx")
    barcode_path = first_existing(sample_dir / "barcodes.tsv.gz", sample_dir / "barcodes.tsv")
    feature_path = first_existing(sample_dir / "features.tsv.gz", sample_dir / "features.tsv")

    # Matrix Market files need mmread; convert to CSC because 10x H5 stores sparse
    # matrices as data/indices/indptr arrays.
    opener = gzip.open if matrix_path.suffix == ".gz" else open
    with opener(matrix_path, "rb") as handle:
        matrix = csc_matrix(mmread(handle))

    barcodes = read_lines(barcode_path)
    features = [line.split("\t") for line in read_lines(feature_path)]

    # 10x matrices are features x barcodes.
    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError(
            f"{sample_dir.name}: matrix shape {matrix.shape} does not match "
            f"{len(features)} features x {len(barcodes)} barcodes"
        )

    feature_ids = [row[0] for row in features]
    feature_names = [row[1] if len(row) > 1 else row[0] for row in features]
    feature_types = [row[2] if len(row) > 2 else "Gene Expression" for row in features]

    # Write the minimal 10x-like H5 structure.
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(output_path, "w") as h5:
        matrix_group = h5.create_group("matrix")
        matrix_group.create_dataset("barcodes", data=np.asarray(barcodes, dtype="S"), compression="gzip")
        matrix_group.create_dataset("data", data=matrix.data.astype("uint32"), compression="gzip")
        matrix_group.create_dataset("indices", data=matrix.indices.astype("uint32"), compression="gzip")
        matrix_group.create_dataset("indptr", data=matrix.indptr.astype("uint32"), compression="gzip")
        matrix_group.create_dataset("shape", data=np.asarray(matrix.shape, dtype="uint64"), compression="gzip")

        features_group = matrix_group.create_group("features")
        features_group.create_dataset("id", data=np.asarray(feature_ids, dtype="S"), compression="gzip")
        features_group.create_dataset("name", data=np.asarray(feature_names, dtype="S"), compression="gzip")
        features_group.create_dataset("feature_type", data=np.asarray(feature_types, dtype="S"), compression="gzip")
        features_group.create_dataset("genome", data=np.asarray([""] * len(features), dtype="S1"), compression="gzip")
        features_group.create_dataset("_all_tag_keys", data=np.asarray([], dtype="S1"), compression="gzip")

    print(f"Wrote {output_path}")


def iter_sample_dirs(input_dir: Path, sample_name: str | None) -> list[Path]:
    """Return selected sample directories for one count type."""
    if sample_name:
        return [input_dir / sample_name]

    if not input_dir.exists():
        return []

    return sorted(path for path in input_dir.iterdir() if path.is_dir())


def convert_count_type(count_type: str, sample_name: str | None) -> int:
    """Convert all selected samples for one count type; return number converted."""
    config = COUNT_CONFIGS[count_type]
    input_dir = config["input_dir"]
    output_dir = config["output_dir"]
    h5_name = config["h5_name"]
    sample_dirs = iter_sample_dirs(input_dir, sample_name)

    if sample_name and not sample_dirs[0].is_dir():
        raise FileNotFoundError(f"{count_type} sample directory not found: {sample_dirs[0]}")

    converted_count = 0
    for sample_dir in sample_dirs:
        output_path = output_dir / sample_dir.name / h5_name
        if is_10x_h5_valid(output_path):
            print(f"Skipping {count_type} {sample_dir.name}: H5 file already exists and is valid.")
            continue

        write_h5(sample_dir, output_path)
        converted_count += 1

    return converted_count


def main() -> None:
    # Convert one sample if SAMPLE_NAME is set; otherwise convert every sample.
    # COUNT_TYPE can be raw, processed, or all.
    sample_name = os.environ.get("SAMPLE_NAME")
    count_type = os.environ.get("COUNT_TYPE", "all").lower()

    if count_type == "all":
        count_types = ["raw", "processed"]
    elif count_type in COUNT_CONFIGS:
        count_types = [count_type]
    else:
        raise ValueError("COUNT_TYPE must be raw, processed, or all")

    total_converted = 0
    for selected_count_type in count_types:
        total_converted += convert_count_type(selected_count_type, sample_name)

    if total_converted == 0:
        print("No H5 files needed conversion.")


if __name__ == "__main__":
    main()
