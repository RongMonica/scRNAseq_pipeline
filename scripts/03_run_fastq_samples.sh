#!/usr/bin/env bash

# Run the FASTQ-to-counts preprocessing pipeline across multiple samples.
#
# This is the batch runner for Path A:
#   FASTQ files -> scripts/00_preprocess_fastq_to_counts.sh
#
# Usage:
#   bash scripts/03_run_fastq_samples.sh
#   SAMPLE_SHEET=samples.csv bash scripts/03_run_fastq_samples.sh
#
# SAMPLE_SHEET can be comma- or tab-separated with columns:
#   sample_name,r1_fastq,r2_fastq

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
DATA_DIR="${DATA_DIR:-$PROJECT_DIR/data}"
FASTQ_DIR="${FASTQ_DIR:-$DATA_DIR}"
SAMPLE_SHEET="${SAMPLE_SHEET:-}"
R1_SUFFIX="${R1_SUFFIX:-_R1.fastq.gz}"
R2_SUFFIX="${R2_SUFFIX:-_R2.fastq.gz}"

found_any="no"

run_sample() {
    local sample_name="$1"
    local r1_fastq="$2"
    local r2_fastq="$3"

    if [[ ! -f "$r1_fastq" ]]; then
        echo "R1 FASTQ not found for $sample_name: $r1_fastq" >&2
        exit 1
    fi

    if [[ ! -f "$r2_fastq" ]]; then
        echo "R2 FASTQ not found for $sample_name: $r2_fastq" >&2
        exit 1
    fi

    echo "Running pipeline for sample: $sample_name"
    SAMPLE_NAME="$sample_name" \
    R1_FASTQ="$r1_fastq" \
    R2_FASTQ="$r2_fastq" \
    bash "$SCRIPT_DIR/00_preprocess_fastq_to_counts.sh"
}

if [[ -n "$SAMPLE_SHEET" ]]; then
    if [[ ! -f "$SAMPLE_SHEET" ]]; then
        echo "Sample sheet not found: $SAMPLE_SHEET" >&2
        exit 1
    fi

    while IFS=$'\t' read -r sample_name r1_fastq r2_fastq extra; do
        if [[ -z "${sample_name:-}" || "$sample_name" == "sample_name" ]]; then
            continue
        fi

        found_any="yes"
        run_sample "$sample_name" "$r1_fastq" "$r2_fastq"
    done < <(tr ',' '\t' < "$SAMPLE_SHEET")
else
    for r1_fastq in "$FASTQ_DIR"/*"$R1_SUFFIX"; do
        if [[ ! -e "$r1_fastq" ]]; then
            continue
        fi

        found_any="yes"
        sample_name="$(basename "$r1_fastq")"
        sample_name="${sample_name%"$R1_SUFFIX"}"
        r2_fastq="$FASTQ_DIR/${sample_name}${R2_SUFFIX}"
        run_sample "$sample_name" "$r1_fastq" "$r2_fastq"
    done
fi

if [[ "$found_any" != "yes" ]]; then
    if [[ -n "$SAMPLE_SHEET" ]]; then
        echo "No samples found in sample sheet: $SAMPLE_SHEET" >&2
    else
        echo "No R1 FASTQ files found in $FASTQ_DIR with suffix $R1_SUFFIX" >&2
    fi
    exit 1
fi
