"""Convert arranged 10x matrix/barcode/feature triplet files into raw H5 files.

Input:
    data/raw/<sample>/matrix.mtx[.gz]
    data/raw/<sample>/barcodes.tsv[.gz]
    data/raw/<sample>/features.tsv[.gz]

Output:
    data/raw_h5/<sample>/raw_feature_bc_matrix.h5

Usage:
    python scripts/02_make_raw_h5.py
    SAMPLE_NAME=<sample> python scripts/02_make_raw_h5.py
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


def write_h5(sample_dir: Path) -> None:
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

    # Write the minimal 10x-like raw_feature_bc_matrix.h5 structure.
    output_path = RAW_H5_DIR / sample_dir.name / "raw_feature_bc_matrix.h5"
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


def main() -> None:
    # Convert one sample if SAMPLE_NAME is set; otherwise convert every sample.
    sample_name = os.environ.get("SAMPLE_NAME")
    sample_dirs = [RAW_DIR / sample_name] if sample_name else sorted(path for path in RAW_DIR.iterdir() if path.is_dir())

    for sample_dir in sample_dirs:
        if not sample_dir.is_dir():
            raise FileNotFoundError(f"Sample directory not found: {sample_dir}")
        write_h5(sample_dir)


if __name__ == "__main__":
    main()
