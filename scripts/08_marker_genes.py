"""Find marker genes, annotate clusters, and run optional downstream analyses.

Input:
    results/clustering/clustered.h5ad

Output:
    results/markers/cluster_markers.tsv
    results/annotation/cluster_annotations.tsv
    results/annotation/annotated.h5ad
    plots/annotation/umap_cell_type.png

Optional:
    results/markers/enrichment.tsv, with --run-enrichment and gseapy
    results/velocity/velocity.h5ad, with --run-velocity and scvelo

Usage:
    python3 scripts/08_marker_genes.py
    python3 scripts/08_marker_genes.py --groupby leiden --run-enrichment
    python3 scripts/08_marker_genes.py --marker-sets data/metadata/cell_type_markers.tsv
"""

import argparse
from collections import defaultdict
from pathlib import Path

import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = REPO_ROOT / "results" / "clustering" / "clustered.h5ad"
MARKERS_DIR = REPO_ROOT / "results" / "markers"
ANNOTATION_DIR = REPO_ROOT / "results" / "annotation"
VELOCITY_DIR = REPO_ROOT / "results" / "velocity"
ANNOTATION_PLOTS_DIR = REPO_ROOT / "plots" / "annotation"

BROAD_AML_MARKERS = {
    "T_cells": ["CD3D", "CD3E", "CD2", "TRAC", "IL7R"],
    "CD8_T_cells": ["CD8A", "CD8B", "GZMK", "NKG7"],
    "NK_cells": ["NKG7", "GNLY", "KLRD1", "PRF1", "GZMB"],
    "B_cells": ["MS4A1", "CD79A", "CD79B", "BANK1"],
    "Plasma_cells": ["MZB1", "JCHAIN", "XBP1", "SDC1"],
    "Monocytes": ["LYZ", "S100A8", "S100A9", "FCN1", "CTSS"],
    "Dendritic_cells": ["FCER1A", "CST3", "CLEC10A", "LILRA4"],
    "HSPC": ["CD34", "SPINK2", "AVP", "PROM1", "GATA2"],
    "Myeloid": ["MPO", "ELANE", "AZU1", "PRTN3", "LYZ"],
    "Erythroid": ["HBB", "HBA1", "HBA2", "GYPA", "ALAS2"],
    "Megakaryocyte": ["PPBP", "PF4", "GP9", "ITGA2B"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find cluster marker genes and assign simple marker-based cell type labels."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=INPUT_PATH,
        help="Clustered AnnData file. Default: results/clustering/clustered.h5ad.",
    )
    parser.add_argument(
        "--groupby",
        default="leiden",
        help="Observation column for cluster groups. Default: leiden.",
    )
    parser.add_argument(
        "--marker-sets",
        type=Path,
        help="Optional TSV/CSV with columns cell_type and gene for annotation markers.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=100,
        help="Number of marker genes saved per cluster. Default: 100.",
    )
    parser.add_argument(
        "--run-enrichment",
        action="store_true",
        help="Run Enrichr gene-set enrichment on top marker genes. Requires gseapy and internet.",
    )
    parser.add_argument(
        "--enrichr-library",
        default="GO_Biological_Process_2023",
        help="Enrichr library used with --run-enrichment. Default: GO_Biological_Process_2023.",
    )
    parser.add_argument(
        "--run-velocity",
        action="store_true",
        help="Run RNA velocity if spliced and unspliced layers are present. Requires scvelo.",
    )
    return parser.parse_args()


def load_marker_sets(path: Path | None) -> dict[str, list[str]]:
    if path is None:
        return BROAD_AML_MARKERS

    sep = "," if path.suffix == ".csv" else "\t"
    table = pd.read_csv(path, sep=sep)
    required = {"cell_type", "gene"}
    if not required.issubset(table.columns):
        raise ValueError(f"Marker set file must contain columns: {', '.join(sorted(required))}")

    marker_sets = defaultdict(list)
    for row in table.itertuples(index=False):
        marker_sets[row.cell_type].append(str(row.gene).upper())
    return dict(marker_sets)


def save_umap(adata: ad.AnnData, color: str, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sc.pl.umap(adata, color=color, frameon=False, show=False)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def find_cluster_markers(adata: ad.AnnData, groupby: str, top_n: int) -> pd.DataFrame:
    if groupby not in adata.obs:
        raise KeyError(f"Groupby column not found in adata.obs: {groupby}")

    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon", use_raw=adata.raw is not None)
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers = markers.sort_values(["group", "scores"], ascending=[True, False])
    return markers.groupby("group", group_keys=False).head(top_n)


def annotate_clusters(
    adata: ad.AnnData,
    groupby: str,
    marker_sets: dict[str, list[str]],
) -> pd.DataFrame:
    score_columns = []
    var_names = {gene.upper(): gene for gene in adata.var_names}
    if adata.raw is not None:
        var_names.update({gene.upper(): gene for gene in adata.raw.var_names})

    for cell_type, genes in marker_sets.items():
        present_genes = [var_names[gene.upper()] for gene in genes if gene.upper() in var_names]
        if len(present_genes) < 2:
            continue

        score_column = f"score_{cell_type}"
        sc.tl.score_genes(adata, gene_list=present_genes, score_name=score_column, use_raw=adata.raw is not None)
        score_columns.append(score_column)

    if not score_columns:
        raise ValueError("No marker sets had at least two genes present in the dataset.")

    cluster_scores = adata.obs.groupby(groupby)[score_columns].mean()
    best_scores = cluster_scores.idxmax(axis=1)
    annotations = pd.DataFrame(
        {
            groupby: cluster_scores.index.astype(str),
            "cell_type": best_scores.str.replace("^score_", "", regex=True).to_numpy(),
            "annotation_score": cluster_scores.max(axis=1).to_numpy(),
        }
    )

    annotation_map = dict(zip(annotations[groupby], annotations["cell_type"], strict=True))
    adata.obs["cell_type"] = adata.obs[groupby].astype(str).map(annotation_map).astype("category")
    return annotations


def run_enrichment(markers: pd.DataFrame, output_path: Path, library: str) -> None:
    try:
        import gseapy as gp
    except ImportError as error:
        raise ImportError("Gene-set enrichment requires gseapy. Install gseapy or omit --run-enrichment.") from error

    rows = []
    for cluster, cluster_markers in markers.groupby("group"):
        genes = cluster_markers["names"].dropna().astype(str).head(100).tolist()
        if len(genes) < 5:
            continue

        result = gp.enrichr(
            gene_list=genes,
            gene_sets=library,
            organism="Human",
            outdir=None,
            no_plot=True,
        )
        table = result.results.copy()
        table.insert(0, "cluster", cluster)
        rows.append(table)

    if rows:
        pd.concat(rows, ignore_index=True).to_csv(output_path, sep="\t", index=False)


def run_velocity(adata: ad.AnnData, output_path: Path) -> None:
    if not {"spliced", "unspliced"}.issubset(adata.layers.keys()):
        print("Skipping RNA velocity: spliced and unspliced layers are not present.")
        return

    try:
        import scvelo as scv
    except ImportError as error:
        raise ImportError("RNA velocity requires scvelo. Install scvelo or omit --run-velocity.") from error

    velocity = adata.copy()
    scv.pp.filter_and_normalize(velocity)
    scv.pp.moments(velocity)
    scv.tl.velocity(velocity)
    scv.tl.velocity_graph(velocity)
    velocity.write_h5ad(output_path)


def main() -> None:
    args = parse_args()

    input_path = args.input.expanduser().resolve()
    marker_sets_path = args.marker_sets.expanduser().resolve() if args.marker_sets else None

    MARKERS_DIR.mkdir(parents=True, exist_ok=True)
    ANNOTATION_DIR.mkdir(parents=True, exist_ok=True)
    VELOCITY_DIR.mkdir(parents=True, exist_ok=True)
    ANNOTATION_PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(input_path)
    markers = find_cluster_markers(adata, args.groupby, args.top_n)
    annotations = annotate_clusters(adata, args.groupby, load_marker_sets(marker_sets_path))

    markers_path = MARKERS_DIR / "cluster_markers.tsv"
    annotations_path = ANNOTATION_DIR / "cluster_annotations.tsv"
    annotated_path = ANNOTATION_DIR / "annotated.h5ad"

    markers.to_csv(markers_path, sep="\t", index=False)
    annotations.to_csv(annotations_path, sep="\t", index=False)
    adata.write_h5ad(annotated_path)

    save_umap(adata, args.groupby, ANNOTATION_PLOTS_DIR / f"umap_{args.groupby}.png")
    save_umap(adata, "cell_type", ANNOTATION_PLOTS_DIR / "umap_cell_type.png")

    if args.run_enrichment:
        run_enrichment(markers, MARKERS_DIR / "enrichment.tsv", args.enrichr_library)

    if args.run_velocity:
        run_velocity(adata, VELOCITY_DIR / "velocity.h5ad")

    print(f"Wrote marker genes to {markers_path}")
    print(f"Wrote cluster annotations to {annotations_path}")
    print(f"Wrote annotated object to {annotated_path}")


if __name__ == "__main__":
    main()
