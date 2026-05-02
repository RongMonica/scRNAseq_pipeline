"""Create per-sample QC-ready AnnData files from 10x-style H5 matrices.

Default input priority:
    1. data/filtered_h5/<sample>/filtered_feature_bc_matrix.h5
       - use directly for QC
    2. data/clean_data/<sample>/<sample>_cellbender.h5
       - use existing CellBender output for QC

With PREFER_CELLBENDER=yes:
    1. data/clean_data/<sample>/<sample>_cellbender.h5
       - use existing CellBender output for QC

Output:
    results/qc/<sample>.h5ad

Usage:
    python scripts/05_qc.py

Prefer existing CellBender output when both filtered and CellBender H5 files
are available:
    PREFER_CELLBENDER=yes python scripts/05_qc.py
"""

from pathlib import Path
import os

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
CLEAN_DATA_DIR = DATA_DIR / "clean_data"
FILTERED_H5_DIR = DATA_DIR / "filtered_h5"
PREPROCESS_DIR = REPO_ROOT / "results" / "preprocess"
OUTPUT_DIR = REPO_ROOT / "results" / "qc"
PREFER_CELLBENDER = os.environ.get("PREFER_CELLBENDER", "no").lower() in {
    "1",
    "true",
    "yes",
    "y",
}

def decode_strings(values) -> list[str]:
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]

def sample_id_from_h5_path(path: Path) -> str:
    # STARsolo paths are usually:
    # results/preprocess/<sample>/Solo.out/Gene/filtered/filtered_feature_bc_matrix.h5
    # For this layout, the sample name is several parents above the H5 file.
    if (
        path.parent.name == "filtered"
        and len(path.parents) > 4
        and path.parents[4] == PREPROCESS_DIR
    ):
        return path.parents[3].name

    # For existing _h5 files that transfered from triplets from 01_arrange_triplet_counts.sh (Path B)
    # data/filtered_h5/<sample>/..., the sample name is simply the parent
    # directory.
    return path.parent.name

def find_filtered_h5_files() -> dict[str, Path]:
    sample_files: dict[str, Path] = {}

    # Search specific H5 folders first, then broader pipeline output folders.
    # setdefault keeps the first match for each sample.
    search_roots = [FILTERED_H5_DIR, PREPROCESS_DIR, DATA_DIR]

    for root in search_roots:
        if not root.exists():
            continue
        for path in sorted(root.rglob("filtered_feature_bc_matrix.h5")):
            sample_files.setdefault(sample_id_from_h5_path(path), path)
    return sample_files

def find_cellbender_h5_files() -> dict[str, Path]:
    sample_files: dict[str, Path] = {}

    if not CLEAN_DATA_DIR.exists():
        return sample_files

    for path in sorted(CLEAN_DATA_DIR.glob("*/*_cellbender.h5")):
        sample_id = path.name.removesuffix("_cellbender.h5")
        sample_files[sample_id] = path
    return sample_files

def get_sample_files() -> dict[str, Path]:
    filtered_files = find_filtered_h5_files()
    cellbender_files = find_cellbender_h5_files()
    sample_files: dict[str, Path] = {}

    if PREFER_CELLBENDER:
        sample_files.update(cellbender_files)
        for sample_id, filtered_path in filtered_files.items():
            sample_files.setdefault(sample_id, filtered_path)
    else:
        sample_files.update(filtered_files)
        for sample_id, cellbender_path in cellbender_files.items():
            sample_files.setdefault(sample_id, cellbender_path)

    return sample_files

def read_10x_h5(h5_path: Path) -> ad.AnnData:
    with h5py.File(h5_path, "r") as handle:
        matrix_group = handle["matrix"]
        matrix = csc_matrix(
            (
                matrix_group["data"][:],
                matrix_group["indices"][:],
                matrix_group["indptr"][:],
            ),
            shape=tuple(matrix_group["shape"][:]),
        ).T.tocsr()

        barcodes = decode_strings(matrix_group["barcodes"][:])

        features_group = matrix_group["features"]
        gene_ids = decode_strings(features_group["id"][:])
        gene_symbols = decode_strings(features_group["name"][:])
        feature_types = decode_strings(features_group["feature_type"][:])

    adata = ad.AnnData(X=matrix)
    adata.obs_names = pd.Index(barcodes, dtype="string")
    adata.var = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "gene_symbol": gene_symbols,
            "feature_type": feature_types,
        }
    )
    adata.var_names = pd.Index(gene_symbols, dtype="string")
    adata.var_names_make_unique()
    return adata

def load_sample(sample_id: str, h5_path: Path) -> ad.AnnData:
    adata = read_10x_h5(h5_path)
    mt_mask = adata.var_names.str.upper().str.startswith("MT-")

    total_counts = np.asarray(adata.X.sum(axis=1)).ravel()
    n_genes_by_counts = np.asarray((adata.X > 0).sum(axis=1)).ravel()

    if mt_mask.any():
        mt_counts = np.asarray(adata[:, mt_mask].X.sum(axis=1)).ravel()
        pct_counts_mt = np.divide(
            mt_counts,
            total_counts,
            out=np.zeros_like(total_counts, dtype=float),
            where=total_counts > 0,
        ) * 100
    else:
        pct_counts_mt = np.zeros(adata.n_obs, dtype=float)

    adata.obs["sample_id"] = sample_id
    adata.obs["total_counts"] = total_counts
    adata.obs["n_genes_by_counts"] = n_genes_by_counts
    adata.obs["pct_counts_mt"] = pct_counts_mt
    adata.var["mt"] = mt_mask.to_numpy()
    adata.uns["source_h5"] = str(h5_path)
    return adata

def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    sample_files = get_sample_files()
    if not sample_files:
        raise FileNotFoundError(
            "No filtered_feature_bc_matrix.h5 or existing CellBender H5 files found"
        )

    for sample_id, h5_path in sample_files.items():
        adata = load_sample(sample_id, h5_path)
        adata.write_h5ad(OUTPUT_DIR / f"{sample_id}.h5ad")
        print(f"Wrote {sample_id}.h5ad from {h5_path.name} with shape {adata.shape}")


if __name__ == "__main__":
    main()
