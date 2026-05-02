#!/usr/bin/env bash

set -euo pipefail

# Path A: Convert one pair of 10x FASTQ files into STARsolo count matrices.
#
# Usage:
#   SAMPLE_NAME=AML_sample01 \
#   R1_FASTQ=/path/to/sample_R1.fastq.gz \
#   R2_FASTQ=/path/to/sample_R2.fastq.gz \
#   bash scripts/00_preprocess_fastq_to_counts.sh

SECONDS=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"

SAMPLE_NAME="${SAMPLE_NAME:-sample01}"
R1_FASTQ="${R1_FASTQ:-$PROJECT_DIR/data/${SAMPLE_NAME}_R1.fastq.gz}"
R2_FASTQ="${R2_FASTQ:-$PROJECT_DIR/data/${SAMPLE_NAME}_R2.fastq.gz}"

REFERENCE_DIR="${REFERENCE_DIR:-$PROJECT_DIR/reference}"
STAR_INDEX_DIR="${STAR_INDEX_DIR:-$REFERENCE_DIR/star_index}"

# GTF_FILE defaults to reference/genes.gtf. Use the same Ensembl GTF as the
# matching STAR index/genome FASTA.
#   human: wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
#   mouse: wget ftp://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz
# Then decompress or symlink the matching file to reference/genes.gtf, or set
# GTF_FILE=/path/to/file.gtf. Do not mix species, genome builds, or releases
# between genome.fa, genes.gtf, and star_index/.
GTF_FILE="${GTF_FILE:-$REFERENCE_DIR/genes.gtf}"

# Barcode whitelist for 10x v3 chemistry. Download with:
#   wget -O reference/3M-february-2018.txt.gz https://github.com/noamteyssier/10x_whitelist_mirror/raw/main/3M-february-2018.txt.gz
WHITELIST="${WHITELIST:-$REFERENCE_DIR/3M-february-2018.txt.gz}"

RESULTS_DIR="${RESULTS_DIR:-$PROJECT_DIR/results}"

# This step turns raw FASTQ files into STARsolo alignment/count outputs.
# Keep those upstream preprocessing results separate from downstream QC,
# clustering, annotation, and marker results.
PREPROCESS_DIR="${PREPROCESS_DIR:-$RESULTS_DIR/preprocess}"
OUTPUT_DIR="${OUTPUT_DIR:-$PREPROCESS_DIR/$SAMPLE_NAME}"

THREADS="${THREADS:-8}"
CHEMISTRY="${CHEMISTRY:-10xv3}"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

require_command() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "Required command not found: $1" >&2
        exit 1
    fi
}

require_file() {
    if [[ ! -f "$1" ]]; then
        echo "Required file not found: $1" >&2
        exit 1
    fi
}

require_dir() {
    if [[ ! -d "$1" ]]; then
        echo "Required directory not found: $1" >&2
        exit 1
    fi
}

case "$CHEMISTRY" in
    10xv2)
        CB_LEN=16
        UMI_LEN=10
        ;;
    10xv3)
        CB_LEN=16
        UMI_LEN=12
        ;;
    *)
        echo "Unsupported CHEMISTRY: $CHEMISTRY" >&2
        echo "Use CHEMISTRY=10xv2 or CHEMISTRY=10xv3." >&2
        exit 1
        ;;
esac

require_command STAR
require_file "$R1_FASTQ"
require_file "$R2_FASTQ"
require_file "$GTF_FILE"
require_file "$WHITELIST"
require_dir "$STAR_INDEX_DIR"

mkdir -p "$OUTPUT_DIR"
cd "$PROJECT_DIR"

log "Sample: $SAMPLE_NAME"
log "R1 FASTQ: $R1_FASTQ"
log "R2 FASTQ: $R2_FASTQ"
log "STAR index: $STAR_INDEX_DIR"
log "Output: $OUTPUT_DIR"

# STARsolo expects cDNA first and barcode/UMI second.
STAR \
    --runThreadN "$THREADS" \
    --genomeDir "$STAR_INDEX_DIR" \
    --readFilesIn "$R2_FASTQ" "$R1_FASTQ" \
    --readFilesCommand zcat \
    --sjdbGTFfile "$GTF_FILE" \
    --outFileNamePrefix "$OUTPUT_DIR/" \
    --outSAMtype BAM SortedByCoordinate \
    --soloType CB_UMI_Simple \
    --soloCBstart 1 \
    --soloCBlen "$CB_LEN" \
    --soloUMIstart 17 \
    --soloUMIlen "$UMI_LEN" \
    --soloCBwhitelist "$WHITELIST" \
    --soloFeatures Gene \
    --soloUMIdedup 1MM_CR \
    --soloCellFilter EmptyDrops_CR

log "Filtered matrix: $OUTPUT_DIR/Solo.out/Gene/filtered/"
log "Raw matrix: $OUTPUT_DIR/Solo.out/Gene/raw/"
log "$((SECONDS / 60)) minutes and $((SECONDS % 60)) seconds elapsed"
