from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
CLEAN_DATA_DIR = DATA_DIR / "clean_data"
OUTPUT_DIR = REPO_ROOT / "results" / "qc"


def decode_strings(values) -> list[str]:
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def get_sample_files() -> dict[str, Path]:
    sample_files: dict[str, Path] = {}

    for path in sorted(CLEAN_DATA_DIR.glob("*/*_cellbender_filtered.h5")):
        sample_id = path.name.removesuffix("_cellbender_filtered.h5")
        sample_files[sample_id] = path

    for path in sorted(CLEAN_DATA_DIR.glob("*/*_cellbender.h5")):
        sample_id = path.name.removesuffix("_cellbender.h5")
        sample_files.setdefault(sample_id, path)

    return sample_files


def read_cellbender_h5(h5_path: Path) -> ad.AnnData:
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
    adata = read_cellbender_h5(h5_path)
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
    adata.uns["cellbender_h5"] = str(h5_path)
    return adata


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    sample_files = get_sample_files()
    if not sample_files:
        raise FileNotFoundError(f"No CellBender output files found in {CLEAN_DATA_DIR}")

    for sample_id, h5_path in sample_files.items():
        adata = load_sample(sample_id, h5_path)
        adata.write_h5ad(OUTPUT_DIR / f"{sample_id}.h5ad")
        print(f"Wrote {sample_id}.h5ad from {h5_path.name} with shape {adata.shape}")


if __name__ == "__main__":
    main()
