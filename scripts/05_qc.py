"""Create per-sample QC-ready AnnData files from 10x-style H5 matrices.

Input:
    User-provided directory containing 10x-style .h5 files.

Output:
    results/qc/<sample>.h5ad
    results/qc/qc_summary.tsv
    plots/qc/<sample>_qc.png

Usage:
    python3 scripts/05_qc.py data/filtered_h5
    python3 scripts/05_qc.py data/clean_data
    python3 scripts/05_qc.py /path/to/h5_directory --output-dir results/qc
"""

import argparse
from pathlib import Path

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "results" / "qc"
PLOTS_DIR = REPO_ROOT / "plots" / "qc"


def decode_strings(values) -> list[str]:
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create per-sample QC-ready AnnData files from 10x-style H5 matrices."
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing .h5 files to process. Search is recursive within this directory only.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Directory for output .h5ad files. Default: results/qc.",
    )
    parser.add_argument(
        "--plots-dir",
        type=Path,
        default=PLOTS_DIR,
        help="Directory for output QC plots. Default: plots/qc.",
    )
    return parser.parse_args()


def sample_id_from_h5_path(path: Path) -> str:
    # STARsolo paths are usually:
    # results/preprocess/<sample>/Solo.out/Gene/filtered/filtered_feature_bc_matrix.h5
    # For this layout, the sample name is several parents above the H5 file.
    if (
        path.name in {"filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5"}
        and path.parent.name in {"filtered", "raw"}
        and len(path.parents) > 3
        and path.parents[2].name == "Solo.out"
    ):
        return path.parents[3].name

    if path.name.endswith("_cellbender.h5"):
        return path.name.removesuffix("_cellbender.h5")

    if path.name in {"filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5"}:
        return path.parent.name

    return path.stem


def find_h5_files(input_dir: Path) -> dict[str, Path]:
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    sample_files: dict[str, Path] = {}
    duplicate_samples: dict[str, list[Path]] = {}

    for path in sorted(input_dir.rglob("*.h5")):
        sample_id = sample_id_from_h5_path(path)
        if sample_id in sample_files:
            duplicate_samples.setdefault(sample_id, [sample_files[sample_id]]).append(path)
            continue
        sample_files[sample_id] = path

    if duplicate_samples:
        details = "\n".join(
            f"{sample_id}: {', '.join(str(path) for path in paths)}"
            for sample_id, paths in sorted(duplicate_samples.items())
        )
        raise ValueError(
            "Multiple .h5 files resolved to the same sample ID. "
            "Point to a narrower directory or rename files:\n"
            f"{details}"
        )

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

    obs = pd.DataFrame(index=pd.Index(barcodes))
    var = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "gene_symbol": gene_symbols,
            "feature_type": feature_types,
        },
        index=pd.Index(gene_symbols),
    )
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
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
    adata.var["mt"] = np.asarray(mt_mask)
    adata.uns["source_h5"] = str(h5_path)
    return adata


def summarize_qc(sample_id: str, adata: ad.AnnData) -> dict[str, object]:
    obs = adata.obs
    return {
        "sample_id": sample_id,
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "median_total_counts": float(np.median(obs["total_counts"])),
        "p05_total_counts": float(np.quantile(obs["total_counts"], 0.05)),
        "p95_total_counts": float(np.quantile(obs["total_counts"], 0.95)),
        "median_n_genes_by_counts": float(np.median(obs["n_genes_by_counts"])),
        "p05_n_genes_by_counts": float(np.quantile(obs["n_genes_by_counts"], 0.05)),
        "p95_n_genes_by_counts": float(np.quantile(obs["n_genes_by_counts"], 0.95)),
        "median_pct_counts_mt": float(np.median(obs["pct_counts_mt"])),
        "pct_cells_mt_gt_20": float((obs["pct_counts_mt"] > 20).mean() * 100),
    }


def write_qc_plot(sample_id: str, adata: ad.AnnData, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    obs = adata.obs
    total_counts = obs["total_counts"]
    n_genes = obs["n_genes_by_counts"]
    pct_mt = obs["pct_counts_mt"]

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    axes[0, 0].hist(total_counts, bins=80, color="#9ecae1", edgecolor="none")
    axes[0, 0].set_title("Count depth")
    axes[0, 0].set_xlabel("Total counts")
    axes[0, 0].set_ylabel("Cells")

    axes[0, 1].hist(n_genes, bins=80, color="#9ecae1", edgecolor="none")
    axes[0, 1].set_title("Detected genes")
    axes[0, 1].set_xlabel("Number of genes")
    axes[0, 1].set_ylabel("Cells")

    ranked_counts = np.sort(total_counts.to_numpy())[::-1]
    axes[1, 0].plot(np.arange(1, len(ranked_counts) + 1), ranked_counts)
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_title("Barcode rank")
    axes[1, 0].set_xlabel("Barcode rank")
    axes[1, 0].set_ylabel("Count depth")

    scatter = axes[1, 1].scatter(
        total_counts,
        n_genes,
        c=pct_mt,
        s=2,
        cmap="viridis",
        alpha=0.5,
    )
    axes[1, 1].set_title("Genes vs count depth")
    axes[1, 1].set_xlabel("Total counts")
    axes[1, 1].set_ylabel("Number of genes")

    cbar = fig.colorbar(scatter, ax=axes[1, 1])
    cbar.set_label("Mitochondrial counts (%)")

    fig.suptitle(sample_id)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    plots_dir = args.plots_dir.expanduser().resolve()

    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    sample_files = find_h5_files(input_dir)
    if not sample_files:
        raise FileNotFoundError(f"No .h5 files found in {input_dir}")

    summary_rows = []
    for sample_id, h5_path in sample_files.items():
        adata = load_sample(sample_id, h5_path)

        adata.write_h5ad(output_dir / f"{sample_id}.h5ad")
        write_qc_plot(sample_id, adata, plots_dir / f"{sample_id}_qc.png")
        summary_rows.append(summarize_qc(sample_id, adata))

        print(f"Wrote {sample_id}.h5ad from {h5_path.name} with shape {adata.shape}")

    summary_path = output_dir / "qc_summary.tsv"
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
    print(f"Wrote QC summary to {summary_path}")


if __name__ == "__main__":
    main()
