#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

if [[ "${1:-}" == "--help" ]]; then
  cat <<'USAGE'
Usage: bash scripts/run_pipeline.sh [--target TARGET]

Builds the dependency-tracked downstream pipeline. Omit --target to build all
analysis targets. Existing artifacts remain at rds/, table/, and figures/.
USAGE
  exit 0
fi

exec Rscript --vanilla code/pipeline/run.R "$@"
