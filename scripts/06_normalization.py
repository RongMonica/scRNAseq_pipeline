"""Filter QC-ready AnnData files and normalize counts per sample.

Input:
    results/qc/<sample>.h5ad
    results/qc/qc_thresholds.tsv, optional

Output:
    results/normalized/<sample>.h5ad

Usage:
    python3 scripts/06_normalization.py
    python3 scripts/06_normalization.py --input-dir results/qc --output-dir results/normalized
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = REPO_ROOT / "results" / "qc"
OUTPUT_DIR = REPO_ROOT / "results" / "normalized"
THRESHOLDS_PATH = INPUT_DIR / "qc_thresholds.tsv"
TARGET_SUM = 10_000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter QC-ready .h5ad files and normalize counts per cell."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=INPUT_DIR,
        help="Directory containing QC .h5ad files. Default: results/qc.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Directory for normalized .h5ad files. Default: results/normalized.",
    )
    parser.add_argument(
        "--thresholds",
        type=Path,
        default=THRESHOLDS_PATH,
        help="QC threshold TSV. Default: results/qc/qc_thresholds.tsv.",
    )
    parser.add_argument(
        "--target-sum",
        type=float,
        default=TARGET_SUM,
        help="Counts per cell after library-size normalization. Default: 10000.",
    )
    return parser.parse_args()


def load_thresholds(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t").set_index("sample_id")


def get_threshold(sample_id: str, thresholds: pd.DataFrame) -> pd.Series | None:
    if thresholds.empty:
        return None
    if sample_id in thresholds.index:
        return thresholds.loc[sample_id]
    if "default" in thresholds.index:
        return thresholds.loc["default"]
    return None


def qc_filter(adata: ad.AnnData, threshold: pd.Series | None) -> ad.AnnData:
    if threshold is None:
        adata.obs["qc_pass"] = True
        return adata.copy()

    obs = adata.obs
    qc_pass = (
        (obs["total_counts"] >= threshold["min_counts"])
        & (obs["total_counts"] <= threshold["max_counts"])
        & (obs["n_genes_by_counts"] >= threshold["min_genes"])
        & (obs["n_genes_by_counts"] <= threshold["max_genes"])
        & (obs["pct_counts_mt"] <= threshold["max_pct_mt"])
    )

    adata.obs["qc_pass"] = qc_pass.to_numpy()
    return adata[qc_pass].copy()


def normalize_log1p(adata: ad.AnnData, target_sum: float) -> ad.AnnData:
    adata.layers["counts"] = adata.X.copy()

    counts_per_cell = np.asarray(adata.X.sum(axis=1)).ravel()
    scale = np.divide(
        target_sum,
        counts_per_cell,
        out=np.zeros_like(counts_per_cell, dtype=float),
        where=counts_per_cell > 0,
    )

    if sparse.issparse(adata.X):
        adata.X = adata.X.multiply(scale[:, None]).tocsr()
        adata.X.data = np.log1p(adata.X.data)
    else:
        adata.X = np.log1p(adata.X * scale[:, None])

    adata.uns["normalization"] = {
        "method": "library_size_normalize_log1p",
        "target_sum": target_sum,
    }
    return adata


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    thresholds_path = args.thresholds.expanduser().resolve()

    output_dir.mkdir(parents=True, exist_ok=True)

    h5ad_files = sorted(input_dir.glob("*.h5ad"))
    if not h5ad_files:
        raise FileNotFoundError(f"No .h5ad files found in {input_dir}")

    thresholds = load_thresholds(thresholds_path)

    for h5ad_path in h5ad_files:
        sample_id = h5ad_path.stem
        adata = ad.read_h5ad(h5ad_path)
        n_cells_before = adata.n_obs

        adata = qc_filter(adata, get_threshold(sample_id, thresholds))
        adata = normalize_log1p(adata, args.target_sum)
        adata.write_h5ad(output_dir / f"{sample_id}.h5ad")

        print(
            f"Wrote {sample_id}.h5ad: "
            f"{n_cells_before} cells before QC, {adata.n_obs} cells after QC"
        )


if __name__ == "__main__":
    main()
