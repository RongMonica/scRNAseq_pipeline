#!/usr/bin/env bash

# Arrange raw matrix/barcode/feature triplet files from a flat input directory
# into per-sample directories expected by downstream preprocessing.
#
# Path B: use this script when starting from existing count triplet files.
#   flat triplet files -> data/raw/<sample>/ -> raw H5 / downstream analysis
#
# Expected input files in SOURCE_DIR:
#   <sample>_raw_matrix.mtx[.gz]
#   <sample>_raw_barcodes.tsv[.gz]
#   <sample>_raw_features.tsv[.gz] or <sample>_raw_genes.tsv[.gz]
#
# Usage:
#   bash scripts/01_arrange_triplet_counts.sh
#   MODE=move bash scripts/01_arrange_triplet_counts.sh
#
# Optional environment variables:
#   SOURCE_DIR=/path/to/flat_raw  # default: data/flat_raw
#   DEST_DIR=/path/to/raw         # default: data/raw
#   MODE=dry-run                  # default; print planned moves
#   MODE=move                     # actually move files

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
SOURCE_DIR="${SOURCE_DIR:-$PROJECT_DIR/data/flat_raw}"
DEST_DIR="${DEST_DIR:-$PROJECT_DIR/data/raw}"
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

declare -A matrices=()
declare -A barcodes=()
declare -A features=()
declare -A samples=()

while IFS= read -r -d '' path; do
    filename="$(basename "$path")"

    case "$filename" in
        *_raw_matrix.mtx.gz)
            sample_name="${filename%_raw_matrix.mtx.gz}"
            matrices["$sample_name"]="$path"
            samples["$sample_name"]=1
            ;;
        *_raw_matrix.mtx)
            sample_name="${filename%_raw_matrix.mtx}"
            matrices["$sample_name"]="$path"
            samples["$sample_name"]=1
            ;;
        *_raw_barcodes.tsv.gz)
            sample_name="${filename%_raw_barcodes.tsv.gz}"
            barcodes["$sample_name"]="$path"
            samples["$sample_name"]=1
            ;;
        *_raw_barcodes.tsv)
            sample_name="${filename%_raw_barcodes.tsv}"
            barcodes["$sample_name"]="$path"
            samples["$sample_name"]=1
            ;;
        *_raw_features.tsv.gz)
            sample_name="${filename%_raw_features.tsv.gz}"
            features["$sample_name"]="$path"
            samples["$sample_name"]=1
            ;;
        *_raw_features.tsv)
            sample_name="${filename%_raw_features.tsv}"
            features["$sample_name"]="$path"
            samples["$sample_name"]=1
            ;;
        *_raw_genes.tsv.gz)
            sample_name="${filename%_raw_genes.tsv.gz}"
            features["$sample_name"]="$path"
            samples["$sample_name"]=1
            ;;
        *_raw_genes.tsv)
            sample_name="${filename%_raw_genes.tsv}"
            features["$sample_name"]="$path"
            samples["$sample_name"]=1
            ;;
        *)
            echo "Skipping unrecognized file: $filename"
            ;;
    esac
done < <(find "$SOURCE_DIR" -maxdepth 1 -type f -print0)

sample_count="${#samples[@]}"
[[ "$sample_count" -gt 0 ]] || fail "No triplet files found in $SOURCE_DIR"

echo "Found $sample_count candidate samples in $SOURCE_DIR"

for sample_name in "${!samples[@]}"; do
    matrix_path="${matrices[$sample_name]:-}"
    barcode_path="${barcodes[$sample_name]:-}"
    feature_path="${features[$sample_name]:-}"

    [[ -n "$matrix_path" ]] || fail "Missing matrix file for $sample_name"
    [[ -n "$barcode_path" ]] || fail "Missing barcode file for $sample_name"
    [[ -n "$feature_path" ]] || fail "Missing feature/gene file for $sample_name"

    sample_dir="$DEST_DIR/$sample_name"
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
    echo "Arrangement complete in $DEST_DIR"
else
    echo "Dry run complete. Re-run with MODE=move to apply changes."
fi
