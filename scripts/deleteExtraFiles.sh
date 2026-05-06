#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
RESULT_DIR="${RESULT_DIR:-$PROJECT_DIR/results/normalized}"


KEEP_SAMPLES=(
  GSM7494260_AML6_DX
  GSM7494261_AML6_REL
  GSM7494262_AML6_REM
  GSM7494295_AML11_DX
  GSM7494296_AML11_REL
  GSM7494297_AML11_REM
  GSM7494257_AML16_DX
  GSM7494258_AML16_REL
  GSM7494259_AML16_REM
)

echo "Deleting files:"
find "${RESULT_DIR}" -maxdepth 1 -type f -name "*.h5ad" \
  ! -name "GSM7494260_AML6_DX.h5ad" \
  ! -name "GSM7494261_AML6_REL.h5ad" \
  ! -name "GSM7494262_AML6_REM.h5ad" \
  ! -name "GSM7494295_AML11_DX.h5ad" \
  ! -name "GSM7494296_AML11_REL.h5ad" \
  ! -name "GSM7494297_AML11_REM.h5ad" \
  ! -name "GSM7494257_AML16_DX.h5ad" \
  ! -name "GSM7494258_AML16_REL.h5ad" \
  ! -name "GSM7494259_AML16_REM.h5ad" \
  -print \
  -delete
