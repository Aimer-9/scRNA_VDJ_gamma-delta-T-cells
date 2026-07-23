#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="$(cd "${ROOT_DIR}/.." && pwd)"

FASTQ_QC_SH="${ROOT_DIR}/src/01_fastq_qc.sh"
GEN_CFG_PY="${ROOT_DIR}/src/02_generate_cellranger_configs.py"
MULTI_SH="${ROOT_DIR}/src/03_cellranger_multi.sh"

pick_first_existing() {
    local path
    for path in "$@"; do
        if [[ -f "${path}" ]]; then
            echo "${path}"
            return 0
        fi
    done
    return 1
}

PARAMS="$(pick_first_existing \
    "${CODE_DIR}/config/params.yaml" \
    "${CODE_DIR}/example/config/params.yaml" \
    "${CODE_DIR}/example/input/params.yaml" || true)"
SAMPLES="$(pick_first_existing \
    "${CODE_DIR}/config/samples.csv" \
    "${CODE_DIR}/example/config/samples.csv" \
    "${CODE_DIR}/example/input/samples.csv" || true)"
TEMPLATE="${ROOT_DIR}/src/multi_template.csv"
OUTDIR=""
QC_OUTDIR=""
THREADS=""

usage() {
    cat <<'USAGE'
Usage:
  bash cellranger.sh <subcommand> [options]

Subcommands:
  dry-run   Validate params.yaml + samples.csv values/paths, then print planned commands
  qc        Run cellranger/src/01_fastq_qc.sh
  multi     Generate Cell Ranger multi configs, then run cellranger multi
  all       Run qc, then multi

Options:
  --params PATH     Path to params.yaml
  --samples PATH    Path to samples.csv
  --template PATH   Path to multi_template.csv
  --outdir DIR      Pipeline output root (overrides params.yaml outdir)
  --qc-outdir DIR   FASTQ QC output directory (default: <outdir>/fastq_qc)
  --threads INT     Threads for FastQC in qc/all (default: cellranger.fastqc_threads in params, or 4)
  -h, --help        Show help

Examples:
  bash code/cellranger/cellranger.sh dry-run
  bash code/cellranger/cellranger.sh qc --params code/config/params.yaml --samples code/config/samples.csv
  bash code/cellranger/cellranger.sh multi --params code/config/params.yaml --samples code/config/samples.csv
  bash code/cellranger/cellranger.sh all --params code/config/params.yaml --samples code/config/samples.csv
USAGE
}

require_file() {
    local path="$1"
    local label="$2"
    if [[ -z "${path}" || ! -f "${path}" ]]; then
        echo "[ERROR] ${label} not found: ${path:-<empty>}" >&2
        exit 1
    fi
}

require_cmd() {
    local cmd="$1"
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: ${cmd}" >&2
        exit 1
    fi
}

read_outdir_from_params() {
    python3 - "$1" <<'PY'
import sys

params = sys.argv[1]
try:
    import yaml
except Exception:
    print("")
    sys.exit(0)

try:
    with open(params, encoding="utf-8") as f:
        d = yaml.safe_load(f) or {}
except Exception:
    print("")
    sys.exit(0)

outdir = d.get("outdir", "")
print(outdir if isinstance(outdir, str) else "")
PY
}

read_fastqc_threads_from_params() {
    python3 - "$1" <<'PY'
import sys

params = sys.argv[1]
try:
    import yaml
except Exception:
    print("")
    sys.exit(0)

try:
    with open(params, encoding="utf-8") as f:
        d = yaml.safe_load(f) or {}
except Exception:
    print("")
    sys.exit(0)

value = (d.get("cellranger") or {}).get("fastqc_threads", "")
if isinstance(value, int):
    print(value)
elif isinstance(value, str) and value.strip().isdigit():
    print(value.strip())
else:
    print("")
PY
}

validate_params_and_samples() {
    python3 - "${PARAMS}" "${SAMPLES}" <<'PY'
import csv
import glob
import os
import shutil
import sys

params_path, samples_path = sys.argv[1], sys.argv[2]

try:
    import yaml
except Exception:
    print("[ERROR] Python package 'yaml' (PyYAML) is required.", file=sys.stderr)
    sys.exit(1)

errors = []
warnings = []

def err(msg):
    errors.append(msg)

def warn(msg):
    warnings.append(msg)

def is_placeholder(value):
    return isinstance(value, str) and value.strip().startswith("/path/to/")

def is_nonempty_str(value):
    return isinstance(value, str) and bool(value.strip())

def check_positive_int(value, key):
    if value is None:
        return
    try:
        if int(value) <= 0:
            err(f"{key} must be > 0")
    except Exception:
        err(f"{key} must be an integer")

# ---- params.yaml ----
if not os.path.isfile(params_path):
    err(f"params.yaml not found: {params_path}")
    cfg = {}
else:
    try:
        with open(params_path, encoding="utf-8") as f:
            cfg = yaml.safe_load(f) or {}
    except Exception as exc:
        err(f"Failed to parse params.yaml: {exc}")
        cfg = {}

if not isinstance(cfg, dict):
    err("Top-level YAML must be a mapping")
    cfg = {}

outdir = cfg.get("outdir")
if not is_nonempty_str(outdir):
    err("Missing or empty key: outdir")
elif is_placeholder(outdir):
    err("outdir still uses placeholder path: /path/to/...")

cellranger = cfg.get("cellranger")
if not isinstance(cellranger, dict):
    err("Missing or invalid key: cellranger")
    cellranger = {}

gex_ref = cellranger.get("reference_gex")
if not is_nonempty_str(gex_ref) or is_placeholder(gex_ref):
    err("Missing or placeholder key: cellranger.reference_gex")
elif not os.path.isdir(gex_ref):
    err(f"cellranger.reference_gex directory not found: {gex_ref}")

check_positive_int(cellranger.get("localcores"), "cellranger.localcores")
check_positive_int(cellranger.get("localmem"), "cellranger.localmem")
check_positive_int(cellranger.get("fastqc_threads"), "cellranger.fastqc_threads")

cr_path = cellranger.get("cellranger_path")
if not is_nonempty_str(cr_path):
    warn("cellranger.cellranger_path not set; runtime will use 'cellranger' from PATH")
else:
    cr_path = cr_path.strip()
    if is_placeholder(cr_path):
        err("cellranger.cellranger_path is placeholder: /path/to/...")
    elif os.path.sep in cr_path:
        if not os.path.isfile(cr_path):
            err(f"cellranger.cellranger_path not found: {cr_path}")
        elif not os.access(cr_path, os.X_OK):
            err(f"cellranger.cellranger_path is not executable: {cr_path}")
    elif shutil.which(cr_path) is None:
        warn(f"cellranger command not found on PATH: {cr_path}")

fastqc_path = cellranger.get("fastqc_path")
if is_nonempty_str(fastqc_path):
    fastqc_path = fastqc_path.strip()
    if is_placeholder(fastqc_path):
        err("cellranger.fastqc_path is placeholder: /path/to/...")
    elif os.path.sep in fastqc_path:
        if not os.path.isfile(fastqc_path):
            err(f"cellranger.fastqc_path not found: {fastqc_path}")
        elif not os.access(fastqc_path, os.X_OK):
            err(f"cellranger.fastqc_path is not executable: {fastqc_path}")
    elif shutil.which(fastqc_path) is None:
        warn(f"fastqc command not found on PATH: {fastqc_path}")

# ---- samples.csv ----
if not os.path.isfile(samples_path):
    err(f"samples.csv not found: {samples_path}")
    rows = []
else:
    try:
        with open(samples_path, newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            fields = reader.fieldnames or []
            rows = list(reader)
    except Exception as exc:
        err(f"Failed to parse samples.csv: {exc}")
        fields = []
        rows = []

required_cols = ["sample_id", "gex_fastq_path", "gex_fastq_prefix"]
missing = [c for c in required_cols if c not in fields]
if missing:
    err("samples.csv missing required columns: " + ", ".join(missing))

optional_libs = [
    ("VDJ-T", "vdj_t_fastq_path", "vdj_t_fastq_prefix"),
    ("VDJ-T-GD", "vdj_t_gd_fastq_path", "vdj_t_gd_fastq_prefix"),
    ("VDJ-B", "vdj_b_fastq_path", "vdj_b_fastq_prefix"),
]
all_libs = [("GEX", "gex_fastq_path", "gex_fastq_prefix")] + [
    t for t in optional_libs if t[1] in fields and t[2] in fields
]

if not rows:
    err("samples.csv has no sample rows")

seen_ids = set()
any_vdj = False
any_vdj_gd = False

for idx, row in enumerate(rows, start=2):
    sid = (row.get("sample_id") or "").strip()
    if not sid:
        err(f"samples.csv line {idx}: sample_id is empty")
        continue
    if sid in seen_ids:
        err(f"samples.csv line {idx}: duplicate sample_id: {sid}")
        continue
    seen_ids.add(sid)

    for lib, path_col, prefix_col in all_libs:
        fastq_dir = (row.get(path_col) or "").strip()
        prefix = (row.get(prefix_col) or "").strip()
        required = (lib == "GEX")
        if required and (not fastq_dir or not prefix):
            err(f"samples.csv line {idx} [{sid}] {lib}: {path_col} and {prefix_col} are required")
            continue
        if (bool(fastq_dir) != bool(prefix)):
            err(f"samples.csv line {idx} [{sid}] {lib}: {path_col} and {prefix_col} must be both set or both empty")
            continue
        if not fastq_dir and not prefix:
            continue

        if lib in ("VDJ-T", "VDJ-T-GD", "VDJ-B"):
            any_vdj = True
        if lib == "VDJ-T-GD":
            any_vdj_gd = True

        if is_placeholder(fastq_dir):
            err(f"samples.csv line {idx} [{sid}] {path_col} is placeholder: {fastq_dir}")
            continue
        if not os.path.isdir(fastq_dir):
            err(f"samples.csv line {idx} [{sid}] {path_col} directory not found: {fastq_dir}")
            continue
        if is_placeholder(prefix):
            err(f"samples.csv line {idx} [{sid}] {prefix_col} is placeholder: {prefix}")
            continue

        files = glob.glob(os.path.join(fastq_dir, f"{prefix}*.fastq.gz"))
        files += glob.glob(os.path.join(fastq_dir, f"{prefix}*.fq.gz"))
        if not files:
            err(f"samples.csv line {idx} [{sid}] {lib}: no FASTQ files matched prefix '{prefix}' in {fastq_dir}")
            continue

        r1 = sum(1 for p in files if ("_R1_" in os.path.basename(p) or "_R1." in os.path.basename(p) or os.path.basename(p).endswith("_1.fastq.gz") or os.path.basename(p).endswith("_1.fq.gz")))
        r2 = sum(1 for p in files if ("_R2_" in os.path.basename(p) or "_R2." in os.path.basename(p) or os.path.basename(p).endswith("_2.fastq.gz") or os.path.basename(p).endswith("_2.fq.gz")))
        if r1 == 0 or r2 == 0:
            err(f"samples.csv line {idx} [{sid}] {lib}: expected both R1 and R2 files")
        elif r1 != r2:
            err(f"samples.csv line {idx} [{sid}] {lib}: R1/R2 file counts differ ({r1} vs {r2})")

if any_vdj:
    vdj_ref = cellranger.get("reference_vdj")
    if not is_nonempty_str(vdj_ref) or is_placeholder(vdj_ref):
        err("Missing or placeholder key: cellranger.reference_vdj (required by VDJ samples)")
    elif not os.path.isdir(vdj_ref):
        err(f"cellranger.reference_vdj directory not found: {vdj_ref}")

if any_vdj_gd:
    primers = cellranger.get("inner_enrichment_primers")
    if not is_nonempty_str(primers) or is_placeholder(primers):
        err("Missing or placeholder key: cellranger.inner_enrichment_primers (required by VDJ-T-GD samples)")
    elif not os.path.isfile(primers):
        err(f"cellranger.inner_enrichment_primers file not found: {primers}")

if warnings:
    print("[WARN] Validation warnings:")
    for w in warnings:
        print(f"  - {w}")

if errors:
    print("[ERROR] Dry-run validation failed:", file=sys.stderr)
    for e in errors:
        print(f"  - {e}", file=sys.stderr)
    sys.exit(1)

print(f"[OK] params.yaml + samples.csv validated:\\n  params:  {params_path}\\n  samples: {samples_path}")
PY
}

run_qc() {
    local validate_only="${1:-false}"
    local cmd=(bash "${FASTQ_QC_SH}"
        --params "${PARAMS}"
        --samples "${SAMPLES}"
        --outdir "${QC_OUTDIR}"
        --threads "${THREADS}")

    if [[ "${validate_only}" == "true" ]]; then
        cmd+=(--validate-only)
    fi

    "${cmd[@]}"
}

run_multi() {
    python3 "${GEN_CFG_PY}" \
        --params "${PARAMS}" \
        --samples "${SAMPLES}" \
        --template "${TEMPLATE}" \
        --outdir "${OUTDIR}"

    bash "${MULTI_SH}" \
        --params "${PARAMS}" \
        --outdir "${OUTDIR}"
}

show_planned_commands() {
    echo "[PLAN] QC command:"
    echo "  bash ${FASTQ_QC_SH} --params ${PARAMS} --samples ${SAMPLES} --outdir ${QC_OUTDIR} --threads ${THREADS}"
    echo "[PLAN] Config generation command:"
    echo "  python3 ${GEN_CFG_PY} --params ${PARAMS} --samples ${SAMPLES} --template ${TEMPLATE} --outdir ${OUTDIR}"
    echo "[PLAN] Cell Ranger multi command:"
    echo "  bash ${MULTI_SH} --params ${PARAMS} --outdir ${OUTDIR}"
}

SUBCOMMAND="${1:-}"
if [[ -z "${SUBCOMMAND}" || "${SUBCOMMAND}" == "-h" || "${SUBCOMMAND}" == "--help" ]]; then
    usage
    exit 0
fi
shift

case "${SUBCOMMAND}" in
    dry-run|qc|multi|all) ;;
    *)
        echo "[ERROR] Unknown subcommand: ${SUBCOMMAND}" >&2
        usage
        exit 1
        ;;
esac

while [[ $# -gt 0 ]]; do
    case "$1" in
        --params)
            PARAMS="$2"
            shift 2
            ;;
        --samples)
            SAMPLES="$2"
            shift 2
            ;;
        --template)
            TEMPLATE="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --qc-outdir)
            QC_OUTDIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            usage
            exit 1
            ;;
    esac
done

require_file "${FASTQ_QC_SH}" "fastq_qc script"
require_file "${GEN_CFG_PY}" "config generation script"
require_file "${MULTI_SH}" "cellranger multi script"

require_file "${PARAMS}" "params.yaml"
require_file "${SAMPLES}" "samples.csv"
require_file "${TEMPLATE}" "multi template"

if [[ -z "${OUTDIR}" ]]; then
    OUTDIR="$(read_outdir_from_params "${PARAMS}")"
fi
if [[ -z "${OUTDIR}" ]]; then
    OUTDIR="output"
fi
if [[ -z "${QC_OUTDIR}" ]]; then
    QC_OUTDIR="${OUTDIR}/fastq_qc"
fi
if [[ -z "${THREADS}" ]]; then
    THREADS="$(read_fastqc_threads_from_params "${PARAMS}")"
fi
if [[ -z "${THREADS}" ]]; then
    THREADS="4"
fi

if ! [[ "${THREADS}" =~ ^[0-9]+$ ]] || [[ "${THREADS}" -lt 1 ]]; then
    echo "[ERROR] --threads must be a positive integer: ${THREADS}" >&2
    exit 1
fi

case "${SUBCOMMAND}" in
    dry-run)
        require_cmd python3
        validate_params_and_samples
        show_planned_commands
        ;;
    qc)
        run_qc false
        ;;
    multi)
        require_cmd Rscript
        run_multi
        ;;
    all)
        require_cmd Rscript
        run_qc false
        run_multi
        ;;
esac
