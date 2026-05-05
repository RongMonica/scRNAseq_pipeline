"""Cluster normalized scRNA-seq AnnData files.

Input:
    results/normalized/<sample>.h5ad

Output:
    results/clustering/clustered.h5ad
    plots/clustering/umap_*.png

Usage:
    python3 scripts/07_clustering.py
    python3 scripts/07_clustering.py --input-dir results/normalized --resolution 0.6
    python3 scripts/07_clustering.py --batch-correction none
"""

import argparse
from pathlib import Path

import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = REPO_ROOT / "results" / "normalized"
OUTPUT_PATH = REPO_ROOT / "results" / "clustering" / "clustered.h5ad"
PLOTS_DIR = REPO_ROOT / "plots" / "clustering"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run PCA, neighbors, UMAP, and Leiden clustering.")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=INPUT_DIR,
        help="Directory containing normalized per-sample .h5ad files. Default: results/normalized.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_PATH,
        help="Output clustered .h5ad file. Default: results/clustering/clustered.h5ad.",
    )
    parser.add_argument(
        "--plots-dir",
        type=Path,
        default=PLOTS_DIR,
        help="Directory for UMAP PNG plots. Default: plots/clustering.",
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=3000,
        help="Number of highly variable genes. Default: 3000.",
    )
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=30,
        help="Number of PCs for neighbor graph. Default: 30.",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=0.6,
        help="Leiden clustering resolution. Default: 0.6.",
    )
    parser.add_argument(
        "--batch-key",
        default="sample_id",
        help="Observation column used as the batch key. Default: sample_id.",
    )
    parser.add_argument(
        "--batch-correction",
        choices=("combat", "none"),
        default="combat",
        help="Batch correction method applied before PCA. Default: combat.",
    )
    return parser.parse_args()


def read_normalized_samples(input_dir: Path) -> ad.AnnData:
    h5ad_files = sorted(
        path
        for path in input_dir.glob("*.h5ad")
        if not path.name.startswith("merged_")
    )
    if not h5ad_files:
        raise FileNotFoundError(f"No normalized .h5ad files found in {input_dir}")

    samples = {}
    for path in h5ad_files:
        sample = ad.read_h5ad(path)
        sample.obs["sample_id"] = path.stem
        samples[path.stem] = sample

    return ad.concat(samples, join="inner", merge="same", index_unique="-")


def save_umap(adata: ad.AnnData, color: str, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sc.pl.umap(adata, color=color, frameon=False, show=False)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def maybe_correct_batch(adata: ad.AnnData, batch_key: str, method: str) -> None:
    if method == "none":
        adata.uns["batch_correction"] = {
            "method": "none",
            "batch_key": batch_key,
        }
        return

    if batch_key not in adata.obs:
        raise KeyError(f"Batch key column not found in adata.obs: {batch_key}")

    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        adata.uns["batch_correction"] = {
            "method": "none",
            "batch_key": batch_key,
            "reason": "fewer_than_two_batches",
        }
        return

    if method == "combat":
        sc.pp.combat(adata, key=batch_key)
        adata.uns["batch_correction"] = {
            "method": "combat",
            "batch_key": batch_key,
            "n_batches": int(n_batches),
        }
        return

    raise ValueError(f"Unsupported batch correction method: {method}")


def main() -> None:
    args = parse_args()

    input_dir = args.input_dir.expanduser().resolve()
    output_path = args.output.expanduser().resolve()
    plots_dir = args.plots_dir.expanduser().resolve()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading normalized samples from {input_dir}...", flush=True)
    adata = read_normalized_samples(input_dir)
    adata.raw = adata

    if args.batch_key not in adata.obs:
        raise KeyError(f"Batch key column not found in adata.obs: {args.batch_key}")

    hvg_batch_key = args.batch_key if adata.obs[args.batch_key].nunique() > 1 else None
    print("Selecting highly variable gens...", flush=True)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=args.n_top_genes,
        flavor="seurat",
        batch_key=hvg_batch_key,
    )

    print("Subsetting to HVGs...", flush=True)
    adata = adata[:, adata.var["highly_variable"]].copy()

    print(f"Running batch correction with method {args.batch_correction}...", flush=True)
    maybe_correct_batch(adata, args.batch_key, args.batch_correction)

    print("Scaling data...", flush=True)
    sc.pp.scale(adata, max_value=10, zero_center=False)

    print("Running PCA...", flush=True)
    sc.tl.pca(adata, svd_solver="arpack")

    print("Building neighbor graph...", flush=True)
    sc.pp.neighbors(adata, n_pcs=args.n_pcs)

    print("Computing UMAP...", flush=True)
    sc.tl.umap(adata)

    print("Running Leiden clustering...", flush=True)
    sc.tl.leiden(adata, resolution=args.resolution, key_added="leiden")

    adata.uns["clustering"] = {
        "n_top_genes": args.n_top_genes,
        "n_pcs": args.n_pcs,
        "resolution": args.resolution,
        "method": "leiden",
        "hvg_batch_key": hvg_batch_key,
    }

    print(f"Writing clustered object to {output_path}...", flush=True)
    adata.write_h5ad(output_path)

    print("Saving UMAP plots...", flush=True)
    save_umap(adata, "leiden", plots_dir / "umap_leiden.png")
    save_umap(adata, "sample_id", plots_dir / "umap_sample_id.png")

    print(f"Wrote clustered object to {output_path}")
    print(f"Wrote UMAP plots to {plots_dir}")


if __name__ == "__main__":
    main()
