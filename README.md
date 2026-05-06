# scRNA-seq Pipeline

This repository contains a single-cell RNA-seq analysis pipeline for AML
samples. It supports two starting points:

- **Path A:** raw FASTQ files that still need alignment and gene counting.
- **Path B:** existing 10x-style matrix/barcode/feature triplet files.

After count generation or conversion, the downstream workflow creates QC-ready
AnnData files, applies sample-level QC thresholds, normalizes counts, clusters
cells, finds marker genes, and assigns simple marker-based cell type labels.

## Repository Layout

```text
data/
  flat/           # downloaded or copied flat triplet files for Path B
  raw/            # arranged raw matrix/barcode/features directories
  raw_h5/         # raw_feature_bc_matrix.h5 files
  filtered_h5/    # filtered_feature_bc_matrix.h5 files
  clean_data/     # CellBender-cleaned H5 files
  metadata/       # sample metadata files
  processed/      # arranged processed triplet files

results/
  preprocess/     # STARsolo outputs from FASTQ processing
  qc/             # QC-ready .h5ad files and qc_summary.tsv
  normalized/     # filtered and normalized per-sample .h5ad files
  clustering/     # merged clustered AnnData object
  markers/        # marker gene and optional enrichment tables
  annotation/     # cluster annotations and annotated AnnData object
  velocity/       # optional RNA velocity output

plots/
  qc/
  clustering/
  annotation/

notebooks/
scripts/
reference/
```

## Requirements

Core Python packages:

```text
anndata
h5py
matplotlib
numpy
pandas
scanpy
scipy
```

Optional tools and packages:

```text
STAR            # FASTQ alignment with STARsolo
FastQC          # optional FASTQ QC
MultiQC         # optional QC report aggregation
cellbender      # optional raw matrix background removal
gseapy          # optional marker gene enrichment
scvelo          # optional RNA velocity
```

## Script Overview

```text
scripts/00_preprocess_fastq_to_counts.sh
  Path A: run STARsolo alignment and UMI gene counting from FASTQ files.

scripts/01_arrange_triplet_counts.sh
  Path B: arrange flat triplet files into data/raw/<sample>/ and
  data/processed/<sample>/.

scripts/02_make_h5.py
  Convert arranged raw and processed triplet files into 10x-style H5 files.

scripts/03_run_fastq_samples.sh
  Batch runner for Path A across many FASTQ samples.

scripts/04_cellbender.py
  Optional CellBender remove-background runner for raw_feature_bc_matrix.h5.

scripts/05_qc.py
  Convert 10x-style H5 files into QC-ready .h5ad files, QC plots, and
  results/qc/qc_summary.tsv.

scripts/05a_write_summary_from_h5ad_files.py
  Rebuild results/qc/qc_summary.tsv from existing QC .h5ad files.

scripts/06_normalization.py
  Apply QC thresholds, preserve raw counts in adata.layers["counts"], and run
  Scanpy normalize_total plus log1p transformation per sample.

scripts/07_clustering.py
  Merge normalized samples, select highly variable genes, optionally correct
  batch effects with ComBat, and run PCA, neighbors, UMAP, and Leiden
  clustering.

scripts/08_marker_genes.py
  Find cluster markers, assign simple marker-based cell type labels, save an
  annotated AnnData object, and optionally run enrichment or RNA velocity.
```

Notebook handoffs:

```text
notebooks/01_qc_review.ipynb
  Review QC summaries and plots, then write results/qc/qc_thresholds.tsv.

notebooks/02_clustering_review.ipynb
  Review clustering and UMAP results.

notebooks/03_marker_annotation.ipynb
  Review marker genes and cell type annotation results.
```

## Quick Start

Use this when you already have filtered 10x-style H5 files:

```bash
python3 scripts/05_qc.py data/filtered_h5
python3 scripts/06_normalization.py
python3 scripts/07_clustering.py
python3 scripts/08_marker_genes.py
```

For raw H5 files that need CellBender first:

```bash
CELLBENDER_EXTRA_ARGS="--expected-cells 5000 --total-droplets-included 20000" \
python3 scripts/04_cellbender.py

python3 scripts/05_qc.py data/clean_data
python3 scripts/06_normalization.py
python3 scripts/07_clustering.py
python3 scripts/08_marker_genes.py
```

## Path A: Starting From FASTQ Files

Use this path when your input data are raw sequencing reads.

```text
FASTQ files
  -> scripts/00_preprocess_fastq_to_counts.sh
  -> STARsolo count matrices
  -> optional raw_feature_bc_matrix.h5 or filtered_feature_bc_matrix.h5
  -> optional scripts/04_cellbender.py for raw H5 files
  -> scripts/05_qc.py
```

For 10x-style data:

```text
R1 = cell barcode + UMI
R2 = cDNA transcript read
```

Run one sample:

```bash
SAMPLE_NAME=AML_sample01 \
R1_FASTQ=/path/to/sample_R1.fastq.gz \
R2_FASTQ=/path/to/sample_R2.fastq.gz \
bash scripts/00_preprocess_fastq_to_counts.sh
```

Common optional settings:

```bash
THREADS=8
CHEMISTRY=10xv3
RUN_FASTQC=yes
RUN_MULTIQC=yes
BUILD_STAR_INDEX=no
```

Expected reference files:

```text
reference/star_index/
reference/genome.fa
reference/genes.gtf
reference/3M-february-2018.txt.gz
```

Main outputs:

```text
results/preprocess/<sample>/Solo.out/Gene/filtered/
results/preprocess/<sample>/Solo.out/Gene/raw/
results/preprocess/<sample>/Aligned.sortedByCoord.out.bam
```

Run many samples with default FASTQ suffixes:

```bash
bash scripts/03_run_fastq_samples.sh
```

Override the FASTQ directory or suffixes:

```bash
FASTQ_DIR=/path/to/fastqs \
R1_SUFFIX=_R1.fastq.gz \
R2_SUFFIX=_R2.fastq.gz \
bash scripts/03_run_fastq_samples.sh
```

Or use a comma- or tab-separated sample sheet:

```bash
SAMPLE_SHEET=samples.csv bash scripts/03_run_fastq_samples.sh
```

```text
sample_name,r1_fastq,r2_fastq
AML_sample01,/path/to/sample01_R1.fastq.gz,/path/to/sample01_R2.fastq.gz
AML_sample02,/path/to/sample02_R1.fastq.gz,/path/to/sample02_R2.fastq.gz
```

## Path B: Starting From Existing Triplet Files

Use this path when you already have 10x-style count files such as matrix,
barcodes, and features or genes files.

```text
flat triplet files
  -> scripts/01_arrange_triplet_counts.sh
  -> data/raw/<sample>/ and data/processed/<sample>/
  -> scripts/02_make_h5.py
  -> data/raw_h5/<sample>/raw_feature_bc_matrix.h5
  -> data/filtered_h5/<sample>/filtered_feature_bc_matrix.h5
```

Expected raw input files in `data/flat/`:

```text
<sample>_raw_matrix.mtx.gz
<sample>_raw_barcodes.tsv.gz
<sample>_raw_features.tsv.gz
```

Expected processed input files in `data/flat/`:

```text
<sample>_processed_matrix.mtx.gz
<sample>_processed_barcodes.tsv.gz
<sample>_processed_feature.tsv.gz
```

Uncompressed files are also supported. The arrange script accepts both singular
and plural feature or gene annotation names, such as `feature.tsv`,
`features.tsv`, `gene.tsv`, and `genes.tsv`.

First, dry-run the arrangement:

```bash
bash scripts/01_arrange_triplet_counts.sh
```

Then move the files:

```bash
MODE=move bash scripts/01_arrange_triplet_counts.sh
```

Optional settings:

```bash
SOURCE_DIR=/path/to/flat
DEST_DIR=/path/to/raw
PROCESSED_DEST_DIR=/path/to/processed
MODE=dry-run
MODE=move
```

After arranging the files, make H5 files:

```bash
python3 scripts/02_make_h5.py
```

Convert one sample:

```bash
SAMPLE_NAME=AML_sample01 python3 scripts/02_make_h5.py
```

Convert only raw or only processed triplets:

```bash
COUNT_TYPE=raw python3 scripts/02_make_h5.py
COUNT_TYPE=processed python3 scripts/02_make_h5.py
```

Input:

```text
data/raw/<sample>/matrix.mtx.gz
data/raw/<sample>/barcodes.tsv.gz
data/raw/<sample>/features.tsv.gz
data/processed/<sample>/matrix.mtx.gz
data/processed/<sample>/barcodes.tsv.gz
data/processed/<sample>/features.tsv.gz
```

Output:

```text
data/raw_h5/<sample>/raw_feature_bc_matrix.h5
data/filtered_h5/<sample>/filtered_feature_bc_matrix.h5
```

## CellBender

Run CellBender when you want to remove background RNA from raw 10x-style H5
matrices.

```bash
CELLBENDER_EXTRA_ARGS="--expected-cells 5000 --total-droplets-included 20000" \
python3 scripts/04_cellbender.py
```

The script searches for `raw_feature_bc_matrix.h5` files under:

```text
data/raw_h5/
results/preprocess/
data/
```

Output:

```text
data/clean_data/<sample>/<sample>_cellbender.h5
```

## QC

Run QC on whichever H5 directory you want to use:

```bash
python3 scripts/05_qc.py data/filtered_h5
python3 scripts/05_qc.py data/clean_data
python3 scripts/05_qc.py /path/to/h5_directory --output-dir results/qc
```

`scripts/05_qc.py` accepts 10x-style `.h5` files from STARsolo, Cell Ranger,
`scripts/02_make_h5.py`, or `scripts/04_cellbender.py`. It searches
recursively inside the input directory and processes every `.h5` file it finds.

Output:

```text
results/qc/<sample>.h5ad
results/qc/qc_summary.tsv
plots/qc/<sample>_qc.png
```

The summary table contains one row per sample with cell count, gene count,
count-depth summaries, detected-gene summaries, and mitochondrial QC summaries.

To rebuild only the summary table from existing QC `.h5ad` files:

```bash
python3 scripts/05a_write_summary_from_h5ad_files.py
```

## QC Review

After running QC, open:

```text
notebooks/01_qc_review.ipynb
```

Use the notebook to review `results/qc/qc_summary.tsv` and
`plots/qc/<sample>_qc.png`, then write:

```text
results/qc/qc_thresholds.tsv
```

The threshold file should include `sample_id` plus these columns:

```text
min_counts
max_counts
min_genes
max_genes
max_pct_mt
```

You can include a `default` row in `sample_id`; it is used for samples without a
sample-specific threshold row.

## Normalization

Run normalization after QC review:

```bash
python3 scripts/06_normalization.py
```

Input:

```text
results/qc/<sample>.h5ad
results/qc/qc_thresholds.tsv
```

Output:

```text
results/normalized/<sample>.h5ad
```

For each sample, the script filters cells using matching thresholds from
`results/qc/qc_thresholds.tsv`. If no threshold file or matching row is found,
all cells are kept. Raw counts are stored in:

```text
adata.layers["counts"]
```

Then the script runs:

```python
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)
```

Useful options:

```bash
python3 scripts/06_normalization.py --target-sum 10000
python3 scripts/06_normalization.py --resume
python3 scripts/06_normalization.py --compression lzf
python3 scripts/06_normalization.py --compression none
```

## Clustering and Batch Correction

Run clustering after normalization:

```bash
python3 scripts/07_clustering.py
```

Input:

```text
results/normalized/
```

Output:

```text
results/clustering/clustered.h5ad
plots/clustering/umap_leiden.png
plots/clustering/umap_sample_id.png
```

By default, samples are merged and `sample_id` is used as the batch key. Highly
variable gene selection is batch-aware when more than one sample is present.
The default batch correction method is Scanpy ComBat.

Default run:

```bash
python3 scripts/07_clustering.py --batch-key sample_id --batch-correction combat
```

Disable batch correction:

```bash
python3 scripts/07_clustering.py --batch-correction none
```

Useful options:

```bash
python3 scripts/07_clustering.py \
  --n-top-genes 3000 \
  --n-pcs 30 \
  --resolution 0.6
```

Batch correction metadata is stored in:

```text
adata.uns["batch_correction"]
```

## Marker Genes and Annotation

Run marker gene detection and cluster annotation after clustering:

```bash
python3 scripts/08_marker_genes.py
```

Input:

```text
results/clustering/clustered.h5ad
```

Output:

```text
results/markers/cluster_markers.tsv
results/annotation/cluster_annotations.tsv
results/annotation/annotated.h5ad
plots/annotation/umap_leiden.png
plots/annotation/umap_cell_type.png
```

The script ranks genes by cluster with Scanpy Wilcoxon testing. It then scores
clusters against broad AML marker sets and writes the best-scoring cell type for
each cluster.

Use a different grouping column:

```bash
python3 scripts/08_marker_genes.py --groupby leiden
```

Save a different number of marker genes per cluster:

```bash
python3 scripts/08_marker_genes.py --top-n 50
```

Provide custom marker sets:

```bash
python3 scripts/08_marker_genes.py --marker-sets data/metadata/cell_type_markers.tsv
```

The custom marker table must be TSV or CSV and contain:

```text
cell_type
gene
```

Optional enrichment analysis:

```bash
python3 scripts/08_marker_genes.py --run-enrichment
python3 scripts/08_marker_genes.py --run-enrichment --enrichr-library GO_Biological_Process_2023
```

This writes:

```text
results/markers/enrichment.tsv
```

Optional RNA velocity:

```bash
python3 scripts/08_marker_genes.py --run-velocity
```

RNA velocity requires `scvelo` and `spliced` and `unspliced` layers in the
clustered AnnData object. If those layers are missing, the script skips velocity.

## Notes

- `scripts/01_arrange_triplet_counts.sh` uses `MODE=dry-run` by default for
  safety.
- Use `MODE=move` only after checking the dry-run output.
- File names should be consistent across matrix, barcode, and feature files so
  the same sample name can be extracted from each file.
- `features.tsv` and `genes.tsv` are treated as feature annotation files.
- Optional enrichment with `gseapy` requires internet access.
