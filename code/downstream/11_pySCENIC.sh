#!/usr/bin/env bash
set -euo pipefail

# pySCENIC workflow for this project.
#
# The default run uses a balanced 10,000-cell subset instead of the full
# Seurat object. This is the practical first-pass size for local pySCENIC:
# large enough to cover the annotated cell types, small enough to avoid the
# common Arboreto/Dask worker failures seen on full 100k-cell matrices.
#
# Pipeline stages:
#   1. Select a reproducible, balanced subset of cells from the annotated
#      Seurat object. Balancing prevents abundant cell types from dominating
#      GRN inference while retaining at least MIN_CELLS_PER_GROUP cells from
#      smaller groups when the requested total permits.
#   2. Export filtered RNA counts as a cells-by-genes CSV for GRN and ctx.
#      Export is performed in chunks to avoid converting the full sparse
#      Seurat matrix to one large dense matrix in R.
#   3. Convert the CSV to an expression loom in Python. This is required
#      because the installed pySCENIC AUCell CLI only accepts loom input.
#   4. Run pySCENIC grn, ctx, and aucell. Each Python step writes a separate
#      log under OUT_DIR/logs.
#   5. Validate that the AUCell result is a real HDF5 loom containing CellID
#      and RegulonsAUC before reporting success.
#
# Prerequisites:
#   - rds/all_seurat_celltype.rds, or another Seurat RDS passed explicitly.
#   - R packages Seurat and Matrix for selection/export.
#   - A dedicated Python environment containing the versions pinned in
#     code/downstream/envs/pyscenic_environment.yml.
#   - A species-matched TF list, cisTarget ranking databases, and motif
#     annotation table. Use the download-ref subcommand to obtain defaults.
#
# Recommended first run:
#   conda activate pyscenic
#   bash code/downstream/11_pySCENIC.sh download-ref \
#     --out-dir pyscenic/ref --species human
#   source pyscenic/ref/pyscenic_reference.env
#   bash code/downstream/11_pySCENIC.sh \
#     --pyscenic-python "$CONDA_PREFIX/bin/python"
#
# Rerun behavior:
#   - The selected-cell list is reused unless --force-selection is supplied.
#   - The CSV export is reused unless --force-export is supplied.
#   - GRN, ctx, and AUCell are rerun whenever this script reaches Python.
#   - The Python runner removes a stale auc_mtx.loom before AUCell and rebuilds
#     expression_matrix.loom if it is absent, older than the CSV, or invalid.
#
# Important outputs:
#   OUT_DIR/selected_cells_<N>_by_<GROUP>.tsv  selected cell IDs and groups
#   OUT_DIR/expression_matrix.csv              GRN/ctx expression input
#   OUT_DIR/expression_matrix.loom             AUCell expression input
#   OUT_DIR/adjacencies.tsv                    inferred TF-target edges
#   OUT_DIR/regulons.csv                       motif-pruned regulons
#   OUT_DIR/auc_mtx.loom                       cell-by-regulon activity output
#   OUT_DIR/logs/*.log                         selection/export/pySCENIC logs
#
# Usage:
#   bash code/downstream/11_pySCENIC.sh [run] [OPTIONS]
#   bash code/downstream/11_pySCENIC.sh download-ref [OPTIONS]
#
# Required for run, as args or env vars:
#   --tf-list PATH                 TF_LIST
#   --ranking-db PATH              RANKING_DBS; repeat for multiple DBs
#   --motif-annotations PATH       MOTIF_ANNOTATIONS
#
# Important run options:
#   --project-dir DIR              PROJECT_DIR [default: /path/to/project]
#   --seurat-rds PATH              SEURAT_RDS [default: rds/all_seurat_celltype.rds]
#   --out-dir DIR                  OUT_DIR [default: pyscenic_10k]
#   --target-cells INT             TARGET_CELLS [default: 10000]
#   --group-by COLUMN              GROUP_BY [default: cell_type]
#   --min-cells-per-group INT      MIN_CELLS_PER_GROUP [default: 200]
#   --loom-chunk-cells INT         LOOM_CHUNK_CELLS [default: 250]
#   --force-selection              Recreate selected-cell list
#   --force-export                 Re-export expression matrix and metadata
#   --pyscenic-python PATH         PYSCENIC_PYTHON [default: python3]
#   --num-workers INT              NUM_WORKERS [default: 16]
#   --grn-num-workers INT          GRN_NUM_WORKERS [default: 1]
#   --grn-method grnboost2|genie3  GRN_METHOD [default: grnboost2]
#   --ctx-mode MODE                CTX_MODE [default: custom_multiprocessing]
#   -h, --help                     Show help

PROJECT_DIR="${PROJECT_DIR:-/path/to/project}"
SEURAT_RDS="${SEURAT_RDS:-rds/all_seurat_celltype.rds}"
OUT_DIR="${OUT_DIR:-pyscenic_10k}"
LOG_DIR="${LOG_DIR:-}"
ASSAY="${ASSAY:-RNA}"
LAYER="${LAYER:-counts}"
TARGET_CELLS="${TARGET_CELLS:-10000}"
GROUP_BY="${GROUP_BY:-cell_type}"
MIN_CELLS_PER_GROUP="${MIN_CELLS_PER_GROUP:-200}"
MIN_COUNTS_PER_GENE="${MIN_COUNTS_PER_GENE:-3}"
MIN_CELLS_PER_GENE="${MIN_CELLS_PER_GENE:-10}"
EXPORT_CHUNK_CELLS="${EXPORT_CHUNK_CELLS:-500}"
LOOM_CHUNK_CELLS="${LOOM_CHUNK_CELLS:-250}"
FORCE_SELECTION="${FORCE_SELECTION:-0}"
FORCE_EXPORT="${FORCE_EXPORT:-0}"
NUM_WORKERS="${NUM_WORKERS:-16}"
GRN_NUM_WORKERS="${GRN_NUM_WORKERS:-1}"
GRN_METHOD="${GRN_METHOD:-grnboost2}"
SEED="${SEED:-1}"
CTX_MODE="${CTX_MODE:-custom_multiprocessing}"
TF_LIST="${TF_LIST:-}"
RANKING_DBS="${RANKING_DBS:-}"
MOTIF_ANNOTATIONS="${MOTIF_ANNOTATIONS:-}"
PYSCENIC_PYTHON="${PYSCENIC_PYTHON:-python3}"
PYSCENIC_CMD="${PYSCENIC_CMD:-pyscenic}"
SELECT_CELLS_R="${SELECT_CELLS_R:-code/downstream/_cache_/pyscenic_select_cells.R}"
EXPORT_R="${EXPORT_R:-code/downstream/_cache_/pyscenic_export.R}"
PYSCENIC_RUN_PY="${PYSCENIC_RUN_PY:-code/downstream/_cache_/pyscenic_run.py}"

RANKING_DB_ARGS=()

# Print the command-line interface using the current environment-derived
# defaults, so the help output reflects values overridden before invocation.
show_help() {
  cat <<EOF
Usage:
  bash code/downstream/11_pySCENIC.sh [run] [OPTIONS]
  bash code/downstream/11_pySCENIC.sh download-ref [OPTIONS]

Required for run, as args or env vars:
  --tf-list PATH                 TF_LIST
  --ranking-db PATH              RANKING_DBS; repeat for multiple DBs
  --motif-annotations PATH       MOTIF_ANNOTATIONS

Important run options:
  --project-dir DIR              PROJECT_DIR [default: $PROJECT_DIR]
  --seurat-rds PATH              SEURAT_RDS [default: $SEURAT_RDS]
  --out-dir DIR                  OUT_DIR [default: $OUT_DIR]
  --target-cells INT             TARGET_CELLS [default: $TARGET_CELLS]
  --group-by COLUMN              GROUP_BY [default: $GROUP_BY]
  --min-cells-per-group INT      MIN_CELLS_PER_GROUP [default: $MIN_CELLS_PER_GROUP]
  --loom-chunk-cells INT         LOOM_CHUNK_CELLS [default: $LOOM_CHUNK_CELLS]
  --force-selection              Recreate selected-cell list
  --force-export                 Re-export expression matrix and metadata
  --pyscenic-python PATH         PYSCENIC_PYTHON [default: $PYSCENIC_PYTHON]
  --num-workers INT              NUM_WORKERS [default: $NUM_WORKERS]
  --grn-num-workers INT          GRN_NUM_WORKERS [default: $GRN_NUM_WORKERS]
  --grn-method grnboost2|genie3  GRN_METHOD [default: $GRN_METHOD]
  --ctx-mode MODE                CTX_MODE [default: $CTX_MODE]
  -h, --help                     Show help
EOF
}

# Require a value after an option such as --project-dir. This prevents a
# missing argument from being interpreted as the next flag.
need_value() {
  local flag="$1"
  local value="${2:-}"
  if [[ -z "$value" ]]; then
    echo "Missing value for $flag" >&2
    exit 1
  fi
}

# Validate required input files early. A zero-byte file is treated as missing
# because interrupted downloads and failed pipeline steps often leave one.
check_file() {
  local path="$1"
  local label="$2"
  if [[ ! -s "$path" ]]; then
    echo "Missing $label: $path" >&2
    exit 1
  fi
}

# Remove URL query parameters before deriving a local download filename.
url_basename() {
  local url="$1"
  url="${url%%\?*}"
  basename "$url"
}

# Download one reference atomically through a temporary file. Existing
# non-empty files are reused unless the caller requests --force.
download_one() {
  local url="$1"
  local output_path="$2"
  local force="$3"
  local tmp_path="${output_path}.tmp"

  if [[ -s "$output_path" && "$force" != "1" ]]; then
    echo "[SKIP] $output_path"
    return 0
  fi

  mkdir -p "$(dirname "$output_path")"
  echo "[GET]  $url"
  if command -v curl >/dev/null 2>&1; then
    curl -L --fail --retry 3 -o "$tmp_path" "$url"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "$tmp_path" "$url"
  else
    echo "Missing curl or wget for reference download." >&2
    exit 1
  fi
  mv "$tmp_path" "$output_path"
}

# Download species-specific pySCENIC resources and write a shell-compatible
# environment file. The ranking databases use mc_v10_clust, not deprecated
# mc9nr databases.
run_download_ref() {
  local ref_out_dir="pyscenic/ref"
  local species="human"
  local genome=""
  local force="0"

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --out-dir) need_value "$1" "${2:-}"; ref_out_dir="$2"; shift 2 ;;
      --species) need_value "$1" "${2:-}"; species="$2"; shift 2 ;;
      --genome) need_value "$1" "${2:-}"; genome="$2"; shift 2 ;;
      --force) force="1"; shift ;;
      -h|--help) show_help; exit 0 ;;
      *) echo "Unknown download-ref argument: $1" >&2; exit 1 ;;
    esac
  done

  if [[ -z "$genome" ]]; then
    case "$species" in
      human) genome="hg38" ;;
      mouse) genome="mm10" ;;
      *) echo "Unsupported species preset: $species" >&2; exit 1 ;;
    esac
  fi

  local tf_url motif_url ranking_urls
  case "${species}:${genome}" in
    human:hg38)
      tf_url="https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt"
      motif_url="https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
      ranking_urls=(
        "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
        "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
      )
      ;;
    mouse:mm10)
      tf_url="https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt"
      motif_url="https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
      ranking_urls=(
        "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
        "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
      )
      ;;
    *) echo "Unsupported species/genome preset: ${species}/${genome}" >&2; exit 1 ;;
  esac

  local tf_path="${ref_out_dir}/$(url_basename "$tf_url")"
  local motif_path="${ref_out_dir}/$(url_basename "$motif_url")"
  local ranking_paths=()
  local url=""

  download_one "$tf_url" "$tf_path" "$force"
  for url in "${ranking_urls[@]}"; do
    local path="${ref_out_dir}/$(url_basename "$url")"
    download_one "$url" "$path" "$force"
    ranking_paths+=("$path")
  done
  download_one "$motif_url" "$motif_path" "$force"

  local env_file="${ref_out_dir}/pyscenic_reference.env"
  {
    printf 'TF_LIST=%q\n' "$tf_path"
    printf 'RANKING_DBS='
    printf '%q ' "${ranking_paths[@]}"
    printf '\n'
    printf 'MOTIF_ANNOTATIONS=%q\n' "$motif_path"
  } > "$env_file"

  echo "Reference env file: $env_file"
}

# Dispatch the optional reference-download subcommand before parsing normal
# analysis options.
if [[ "${1:-run}" == "download-ref" ]]; then
  shift
  run_download_ref "$@"
  exit 0
elif [[ "${1:-}" == "run" ]]; then
  shift
fi

# Parse run options. Every setting also has an uppercase environment-variable
# equivalent, which is useful for scheduled jobs and reusable shell profiles.
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-dir) need_value "$1" "${2:-}"; PROJECT_DIR="$2"; shift 2 ;;
    --seurat-rds) need_value "$1" "${2:-}"; SEURAT_RDS="$2"; shift 2 ;;
    --out-dir) need_value "$1" "${2:-}"; OUT_DIR="$2"; shift 2 ;;
    --log-dir) need_value "$1" "${2:-}"; LOG_DIR="$2"; shift 2 ;;
    --assay) need_value "$1" "${2:-}"; ASSAY="$2"; shift 2 ;;
    --layer) need_value "$1" "${2:-}"; LAYER="$2"; shift 2 ;;
    --target-cells) need_value "$1" "${2:-}"; TARGET_CELLS="$2"; shift 2 ;;
    --group-by) need_value "$1" "${2:-}"; GROUP_BY="$2"; shift 2 ;;
    --min-cells-per-group) need_value "$1" "${2:-}"; MIN_CELLS_PER_GROUP="$2"; shift 2 ;;
    --min-counts-per-gene) need_value "$1" "${2:-}"; MIN_COUNTS_PER_GENE="$2"; shift 2 ;;
    --min-cells-per-gene) need_value "$1" "${2:-}"; MIN_CELLS_PER_GENE="$2"; shift 2 ;;
    --export-chunk-cells) need_value "$1" "${2:-}"; EXPORT_CHUNK_CELLS="$2"; shift 2 ;;
    --loom-chunk-cells) need_value "$1" "${2:-}"; LOOM_CHUNK_CELLS="$2"; shift 2 ;;
    --num-workers) need_value "$1" "${2:-}"; NUM_WORKERS="$2"; shift 2 ;;
    --grn-num-workers) need_value "$1" "${2:-}"; GRN_NUM_WORKERS="$2"; shift 2 ;;
    --grn-method) need_value "$1" "${2:-}"; GRN_METHOD="$2"; shift 2 ;;
    --seed) need_value "$1" "${2:-}"; SEED="$2"; shift 2 ;;
    --ctx-mode) need_value "$1" "${2:-}"; CTX_MODE="$2"; shift 2 ;;
    --pyscenic-python) need_value "$1" "${2:-}"; PYSCENIC_PYTHON="$2"; shift 2 ;;
    --pyscenic-command) need_value "$1" "${2:-}"; PYSCENIC_CMD="$2"; shift 2 ;;
    --tf-list) need_value "$1" "${2:-}"; TF_LIST="$2"; shift 2 ;;
    --ranking-db) need_value "$1" "${2:-}"; RANKING_DB_ARGS+=("$2"); shift 2 ;;
    --ranking-dbs) need_value "$1" "${2:-}"; RANKING_DBS="$2"; shift 2 ;;
    --motif-annotations) need_value "$1" "${2:-}"; MOTIF_ANNOTATIONS="$2"; shift 2 ;;
    --force-selection) FORCE_SELECTION="1"; shift ;;
    --force-export)
      if [[ "${2:-}" =~ ^[01]$ ]]; then
        FORCE_EXPORT="$2"
        shift 2
      else
        FORCE_EXPORT="1"
        shift
      fi
      ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; show_help >&2; exit 1 ;;
  esac
done

# Normalize ranking databases supplied either as a whitespace-separated
# RANKING_DBS value or as repeated --ranking-db arguments.
RANKING_DB_ARRAY=()
if [[ -n "$RANKING_DBS" ]]; then
  read -r -a RANKING_DB_ARRAY <<< "$RANKING_DBS"
fi
if [[ "${#RANKING_DB_ARGS[@]}" -gt 0 ]]; then
  RANKING_DB_ARRAY+=("${RANKING_DB_ARGS[@]}")
fi

# The three reference resource classes are mandatory and must all describe
# the same species/gene-symbol convention as the expression matrix.
if [[ -z "$TF_LIST" ]]; then
  echo "Set --tf-list or TF_LIST." >&2
  exit 1
fi
if [[ "${#RANKING_DB_ARRAY[@]}" -eq 0 ]]; then
  echo "Set --ranking-db/--ranking-dbs or RANKING_DBS." >&2
  exit 1
fi
if [[ -z "$MOTIF_ANNOTATIONS" ]]; then
  echo "Set --motif-annotations or MOTIF_ANNOTATIONS." >&2
  exit 1
fi

# All relative input and output paths below are interpreted from PROJECT_DIR.
cd "$PROJECT_DIR"
mkdir -p "$OUT_DIR"
if [[ -z "$LOG_DIR" ]]; then
  LOG_DIR="$OUT_DIR/logs"
fi
mkdir -p "$LOG_DIR"

# Fail before expensive computation if any project input, helper, reference,
# or Python executable is unavailable.
check_file "$SEURAT_RDS" "Seurat RDS"
check_file "$TF_LIST" "TF list"
check_file "$MOTIF_ANNOTATIONS" "motif annotations"
check_file "$SELECT_CELLS_R" "cell-selection R script"
check_file "$EXPORT_R" "pySCENIC export R script"
check_file "$PYSCENIC_RUN_PY" "pySCENIC Python runner"
for db in "${RANKING_DB_ARRAY[@]}"; do
  check_file "$db" "ranking database"
done

if ! command -v "$PYSCENIC_PYTHON" >/dev/null 2>&1 && [[ ! -x "$PYSCENIC_PYTHON" ]]; then
  echo "Missing pySCENIC Python executable: $PYSCENIC_PYTHON" >&2
  exit 1
fi

# Define all pipeline products in one place. The selected-cell filename records
# the requested sample size and balancing column for easier provenance.
SELECTED_CELLS="$OUT_DIR/selected_cells_${TARGET_CELLS}_by_${GROUP_BY}.tsv"
EXPR_CSV="$OUT_DIR/expression_matrix.csv"
EXPR_LOOM="$OUT_DIR/expression_matrix.loom"
METADATA_TSV="$OUT_DIR/cell_metadata.tsv"
ADJ_TSV="$OUT_DIR/adjacencies.tsv"
REGULONS_CSV="$OUT_DIR/regulons.csv"
AUC_LOOM="$OUT_DIR/auc_mtx.loom"
EXPORT_LOG="$LOG_DIR/pyscenic_export.log"
SELECT_LOG="$LOG_DIR/pyscenic_select_cells.log"

# Stage 1: select cells reproducibly. The R helper allocates a minimum number
# per group first, then distributes remaining capacity among larger groups.
if [[ "$FORCE_SELECTION" == "1" || ! -s "$SELECTED_CELLS" ]]; then
  echo "Selecting balanced pySCENIC cell subset"
  echo "[LOG] $SELECT_LOG"
  Rscript --vanilla "$SELECT_CELLS_R" \
    "$SEURAT_RDS" \
    "$SELECTED_CELLS" \
    "$GROUP_BY" \
    "$TARGET_CELLS" \
    "$MIN_CELLS_PER_GROUP" \
    "$SEED" > "$SELECT_LOG" 2>&1
else
  echo "[SKIP] Existing selected-cell list: $SELECTED_CELLS"
fi

# Stage 2: export selected cells from the requested Seurat assay/layer.
# The literal zero disables additional random downsampling in the exporter,
# because the balanced selected-cell list already defines the exact cohort.
if [[ "$FORCE_EXPORT" == "1" || ! -s "$EXPR_CSV" || ! -s "$METADATA_TSV" ]]; then
  echo "Exporting expression matrix for selected cells"
  echo "[LOG] $EXPORT_LOG"
  Rscript --vanilla "$EXPORT_R" \
    "$SEURAT_RDS" \
    "$EXPR_CSV" \
    "$METADATA_TSV" \
    "$ASSAY" \
    "$LAYER" \
    "$MIN_COUNTS_PER_GENE" \
    "$MIN_CELLS_PER_GENE" \
    0 \
    "$EXPORT_CHUNK_CELLS" \
    "$SELECTED_CELLS" > "$EXPORT_LOG" 2>&1
else
  echo "[SKIP] Existing expression export: $EXPR_CSV"
fi

# Convert ranking database paths to repeated Python CLI arguments.
PYSCENIC_RANKING_ARGS=()
for db in "${RANKING_DB_ARRAY[@]}"; do
  PYSCENIC_RANKING_ARGS+=(--ranking-db "$db")
done

# Stages 3-5 are handled by the Python runner:
#   - GRN and ctx consume expression_matrix.csv.
#   - CSV is converted to expression_matrix.loom in LOOM_CHUNK_CELLS chunks.
#   - AUCell consumes the expression loom and creates auc_mtx.loom.
#   - All major outputs are validated before this command returns success.
"$PYSCENIC_PYTHON" "$PYSCENIC_RUN_PY" \
  --expr-csv "$EXPR_CSV" \
  --expr-loom "$EXPR_LOOM" \
  --loom-chunk-cells "$LOOM_CHUNK_CELLS" \
  --tf-list "$TF_LIST" \
  --motif-annotations "$MOTIF_ANNOTATIONS" \
  --adj-tsv "$ADJ_TSV" \
  --regulons-csv "$REGULONS_CSV" \
  --auc-loom "$AUC_LOOM" \
  --num-workers "$NUM_WORKERS" \
  --grn-num-workers "$GRN_NUM_WORKERS" \
  --grn-method "$GRN_METHOD" \
  --seed "$SEED" \
  --ctx-mode "$CTX_MODE" \
  --pyscenic-command "$PYSCENIC_CMD" \
  --log-dir "$LOG_DIR" \
  "${PYSCENIC_RANKING_ARGS[@]}"

cat <<EOF
pySCENIC complete.

Selected cells:    $SELECTED_CELLS
Expression matrix: $EXPR_CSV
Expression loom:   $EXPR_LOOM
Cell metadata:     $METADATA_TSV
Adjacencies:       $ADJ_TSV
Regulons:          $REGULONS_CSV
AUCell loom:       $AUC_LOOM
Logs:              $LOG_DIR
EOF
