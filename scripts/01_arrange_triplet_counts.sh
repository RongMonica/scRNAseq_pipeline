#!/usr/bin/env bash

# Arrange raw and processed matrix/barcode/feature triplet files from a flat
# input directory into per-sample directories.
#
# Path B: use this script when starting from existing count triplet files.
#   flat raw triplet files -> data/raw/<sample>/ -> raw H5 / downstream analysis
#   flat processed triplet files -> data/processed/<sample>/
#
# Expected input files in SOURCE_DIR:
#   <sample>_raw_matrix.mtx[.gz]
#   <sample>_raw_barcodes.tsv[.gz]
#   <sample>_raw_features.tsv[.gz] or <sample>_raw_genes.tsv[.gz]
#   <sample>_processed_matrix.mtx[.gz]
#   <sample>_processed_barcodes.tsv[.gz]
#   <sample>_processed_feature.tsv[.gz], <sample>_processed_features.tsv[.gz],
#       <sample>_processed_gene.tsv[.gz], or <sample>_processed_genes.tsv[.gz]
#
# Usage:
#   bash scripts/01_arrange_triplet_counts.sh
#   MODE=move bash scripts/01_arrange_triplet_counts.sh
#
# Optional environment variables:
#   SOURCE_DIR=/path/to/flat      # default: data/flat
#   DEST_DIR=/path/to/raw                 # default: data/raw
#   PROCESSED_DEST_DIR=/path/to/processed # default: data/processed
#   MODE=dry-run                  # default; print planned moves
#   MODE=move                     # actually move files

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
SOURCE_DIR="${SOURCE_DIR:-$PROJECT_DIR/data/flat}"
DEST_DIR="${DEST_DIR:-$PROJECT_DIR/data/raw}"
PROCESSED_DEST_DIR="${PROCESSED_DEST_DIR:-$PROJECT_DIR/data/processed}"
MODE="${MODE:-dry-run}"

fail() {
    echo "$*" >&2
    exit 1
}

move_or_print() {
    local source_path="$1"
    local dest_path="$2"

    if [[ "$MODE" == "move" ]]; then
        mv "$source_path" "$dest_path"
        echo "Moved $(basename "$source_path") -> $dest_path"
    else
        echo "Would move $(basename "$source_path") -> $dest_path"
    fi
}

[[ "$MODE" == "dry-run" || "$MODE" == "move" ]] || fail "MODE must be dry-run or move"
[[ -d "$SOURCE_DIR" ]] || fail "Source directory not found: $SOURCE_DIR"
mkdir -p "$DEST_DIR"
mkdir -p "$PROCESSED_DEST_DIR"

declare -A matrices=()
declare -A barcodes=()
declare -A features=()
declare -A count_sets=()
declare -A count_types=()

register_triplet_file() {
    local count_type="$1"
    local sample_name="$2"
    local file_kind="$3"
    local path="$4"
    local key="$count_type|$sample_name"

    case "$file_kind" in
        matrix)
            matrices["$key"]="$path"
            ;;
        barcode)
            barcodes["$key"]="$path"
            ;;
        feature)
            features["$key"]="$path"
            ;;
        *)
            fail "Unknown triplet file kind: $file_kind"
            ;;
    esac

    count_sets["$key"]=1
    count_types["$count_type"]=1
}

while IFS= read -r -d '' path; do
    filename="$(basename "$path")"

    case "$filename" in
        *_raw_matrix.mtx.gz)
            sample_name="${filename%_raw_matrix.mtx.gz}"
            register_triplet_file "raw" "$sample_name" "matrix" "$path"
            ;;
        *_raw_matrix.mtx)
            sample_name="${filename%_raw_matrix.mtx}"
            register_triplet_file "raw" "$sample_name" "matrix" "$path"
            ;;
        *_raw_barcodes.tsv.gz)
            sample_name="${filename%_raw_barcodes.tsv.gz}"
            register_triplet_file "raw" "$sample_name" "barcode" "$path"
            ;;
        *_raw_barcodes.tsv)
            sample_name="${filename%_raw_barcodes.tsv}"
            register_triplet_file "raw" "$sample_name" "barcode" "$path"
            ;;
        *_raw_feature.tsv.gz)
            sample_name="${filename%_raw_feature.tsv.gz}"
            register_triplet_file "raw" "$sample_name" "feature" "$path"
            ;;
        *_raw_feature.tsv)
            sample_name="${filename%_raw_feature.tsv}"
            register_triplet_file "raw" "$sample_name" "feature" "$path"
            ;;
        *_raw_features.tsv.gz)
            sample_name="${filename%_raw_features.tsv.gz}"
            register_triplet_file "raw" "$sample_name" "feature" "$path"
            ;;
        *_raw_features.tsv)
            sample_name="${filename%_raw_features.tsv}"
            register_triplet_file "raw" "$sample_name" "feature" "$path"
            ;;
        *_raw_gene.tsv.gz)
            sample_name="${filename%_raw_gene.tsv.gz}"
            register_triplet_file "raw" "$sample_name" "feature" "$path"
            ;;
        *_raw_gene.tsv)
            sample_name="${filename%_raw_gene.tsv}"
            register_triplet_file "raw" "$sample_name" "feature" "$path"
            ;;
        *_raw_genes.tsv.gz)
            sample_name="${filename%_raw_genes.tsv.gz}"
            register_triplet_file "raw" "$sample_name" "feature" "$path"
            ;;
        *_raw_genes.tsv)
            sample_name="${filename%_raw_genes.tsv}"
            register_triplet_file "raw" "$sample_name" "feature" "$path"
            ;;
        *_processed_matrix.mtx.gz)
            sample_name="${filename%_processed_matrix.mtx.gz}"
            register_triplet_file "processed" "$sample_name" "matrix" "$path"
            ;;
        *_processed_matrix.mtx)
            sample_name="${filename%_processed_matrix.mtx}"
            register_triplet_file "processed" "$sample_name" "matrix" "$path"
            ;;
        *_processed_barcodes.tsv.gz)
            sample_name="${filename%_processed_barcodes.tsv.gz}"
            register_triplet_file "processed" "$sample_name" "barcode" "$path"
            ;;
        *_processed_barcodes.tsv)
            sample_name="${filename%_processed_barcodes.tsv}"
            register_triplet_file "processed" "$sample_name" "barcode" "$path"
            ;;
        *_processed_feature.tsv.gz)
            sample_name="${filename%_processed_feature.tsv.gz}"
            register_triplet_file "processed" "$sample_name" "feature" "$path"
            ;;
        *_processed_feature.tsv)
            sample_name="${filename%_processed_feature.tsv}"
            register_triplet_file "processed" "$sample_name" "feature" "$path"
            ;;
        *_processed_features.tsv.gz)
            sample_name="${filename%_processed_features.tsv.gz}"
            register_triplet_file "processed" "$sample_name" "feature" "$path"
            ;;
        *_processed_features.tsv)
            sample_name="${filename%_processed_features.tsv}"
            register_triplet_file "processed" "$sample_name" "feature" "$path"
            ;;
        *_processed_gene.tsv.gz)
            sample_name="${filename%_processed_gene.tsv.gz}"
            register_triplet_file "processed" "$sample_name" "feature" "$path"
            ;;
        *_processed_gene.tsv)
            sample_name="${filename%_processed_gene.tsv}"
            register_triplet_file "processed" "$sample_name" "feature" "$path"
            ;;
        *_processed_genes.tsv.gz)
            sample_name="${filename%_processed_genes.tsv.gz}"
            register_triplet_file "processed" "$sample_name" "feature" "$path"
            ;;
        *_processed_genes.tsv)
            sample_name="${filename%_processed_genes.tsv}"
            register_triplet_file "processed" "$sample_name" "feature" "$path"
            ;;
        *)
            echo "Skipping unrecognized file: $filename"
            ;;
    esac
done < <(find "$SOURCE_DIR" -maxdepth 1 -type f -print0)

count_set_count="${#count_sets[@]}"
[[ "$count_set_count" -gt 0 ]] || fail "No raw or processed triplet files found in $SOURCE_DIR"

echo "Found $count_set_count candidate count sets in $SOURCE_DIR"

for key in "${!count_sets[@]}"; do
    count_type="${key%%|*}"
    sample_name="${key#*|}"
    matrix_path="${matrices[$key]:-}"
    barcode_path="${barcodes[$key]:-}"
    feature_path="${features[$key]:-}"

    [[ -n "$matrix_path" ]] || fail "Missing $count_type matrix file for $sample_name"
    [[ -n "$barcode_path" ]] || fail "Missing $count_type barcode file for $sample_name"
    [[ -n "$feature_path" ]] || fail "Missing $count_type feature/gene file for $sample_name"

    if [[ "$count_type" == "processed" ]]; then
        sample_dir="$PROCESSED_DEST_DIR/$sample_name"
    else
        sample_dir="$DEST_DIR/$sample_name"
    fi
    matrix_dest="$sample_dir/matrix.mtx${matrix_path##*.mtx}"
    barcode_dest="$sample_dir/barcodes.tsv${barcode_path##*.tsv}"
    feature_dest="$sample_dir/features.tsv${feature_path##*.tsv}"

    if [[ "$MODE" == "move" ]]; then
        mkdir -p "$sample_dir"
    fi

    [[ ! -e "$matrix_dest" ]] || fail "Destination already exists: $matrix_dest"
    [[ ! -e "$barcode_dest" ]] || fail "Destination already exists: $barcode_dest"
    [[ ! -e "$feature_dest" ]] || fail "Destination already exists: $feature_dest"

    move_or_print "$matrix_path" "$matrix_dest"
    move_or_print "$barcode_path" "$barcode_dest"
    move_or_print "$feature_path" "$feature_dest"
done

if [[ "$MODE" == "move" ]]; then
    if [[ -n "${count_types[raw]:-}" ]]; then
        echo "Raw arrangement complete in $DEST_DIR"
    fi
    if [[ -n "${count_types[processed]:-}" ]]; then
        echo "Processed arrangement complete in $PROCESSED_DEST_DIR"
    fi
else
    echo "Dry run complete. Re-run with MODE=move to apply changes."
fi
