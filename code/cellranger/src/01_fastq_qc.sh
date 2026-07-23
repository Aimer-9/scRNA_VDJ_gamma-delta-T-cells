#!/bin/bash
# ============================================================================
# 01_fastq_qc.sh
# Validate FASTQ inputs declared in example/input/samples.csv and run FastQC/MultiQC
# before Cell Ranger.
#
# Usage:
#   bash cellranger/src/01_fastq_qc.sh [OPTIONS]
#
# Options:
#   --sample PATH      Path to samples.csv [default: example/input/samples.csv]
#   --params PATH
#                      Path to unified params.yaml
#                      [default: example/input/params.yaml]
#   --outdir DIR       Output directory for QC reports
#                      [default: example/output/fastq_qc]
#   --threads INT      Threads passed to FastQC
#                      [default: cellranger.fastqc_threads in params YAML, or 4]
#   --validate-only    Validate inputs only; skip FastQC
#                      and MultiQC
#   --dry-run          Force defaults, skip validation, echo planned actions
#   -h, --help         Show this help
#
# MultiQC is run inside the conda environment named 'multiqc'.
# Create it with:
#   conda create -n multiqc -c bioconda -c conda-forge multiqc
# ============================================================================

set -euo pipefail

# ---- Defaults ---------------------------------------------------------------
DEFAULT_SAMPLES_CSV="example/input/samples.csv"
DEFAULT_PARAMS_YAML="example/input/params.yaml"
DEFAULT_OUTDIR="example/output/fastq_qc"
DEFAULT_FASTQC="fastqc"
DEFAULT_THREADS=4

SAMPLES_CSV="${DEFAULT_SAMPLES_CSV}"
PARAMS_YAML="${DEFAULT_PARAMS_YAML}"
OUTDIR="${DEFAULT_OUTDIR}"
FASTQC_ARG=""
FASTQC="${DEFAULT_FASTQC}"
THREADS=""
VALIDATE_ONLY=false
TEST_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample|--samples) SAMPLES_CSV="$2";   shift 2 ;;
        --params) PARAMS_YAML="$2"; shift 2 ;;
        --outdir)        OUTDIR="$2";        shift 2 ;;
        --fastqc-path)   FASTQC_ARG="$2";    shift 2 ;;
        --threads)       THREADS="$2";       shift 2 ;;
        --validate-only) VALIDATE_ONLY=true; shift   ;;
        --dry-run)       TEST_RUN=true;      shift   ;;
        -h|--help)
            sed -n '/^# Usage/,/^# ====/p' "$0" | sed 's/^# \?//'
            exit 0 ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            exit 1
            ;;
    esac
done

_yaml_get() {
    python3 - "$PARAMS_YAML" "$1" "$2" <<'PY'
import sys

params_yaml, dotted_key, default = sys.argv[1:]
try:
    import yaml
except Exception:
    print(default)
    sys.exit(0)

try:
    with open(params_yaml, encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
except Exception:
    print(default)
    sys.exit(0)

value = data
for key in dotted_key.split("."):
    if not isinstance(value, dict) or key not in value or value[key] is None:
        print(default)
        sys.exit(0)
    value = value[key]

print(value)
PY
}

if [[ -n "${FASTQC_ARG}" ]]; then
    FASTQC="${FASTQC_ARG}"
else
    FASTQC=$(_yaml_get "cellranger.fastqc_path" "${DEFAULT_FASTQC}")
fi
if [[ -z "${THREADS}" ]]; then
    THREADS=$(_yaml_get "cellranger.fastqc_threads" "${DEFAULT_THREADS}")
fi

if [[ "${TEST_RUN}" == true ]]; then
    SAMPLES_CSV="${DEFAULT_SAMPLES_CSV}"
    PARAMS_YAML="${DEFAULT_PARAMS_YAML}"
    OUTDIR="${DEFAULT_OUTDIR}"
    FASTQC_ARG=""
    FASTQC=$(_yaml_get "cellranger.fastqc_path" "${DEFAULT_FASTQC}")
    THREADS=$(_yaml_get "cellranger.fastqc_threads" "${DEFAULT_THREADS}")
    VALIDATE_ONLY=false
    echo "[TEST]  --dry-run: using defaults, skipping validation."
    echo "[TEST]  SAMPLE_CSV  = ${SAMPLES_CSV}"
    echo "[TEST]  PARAMS_YAML = ${PARAMS_YAML}"
    echo "[TEST]  OUTDIR      = ${OUTDIR}"
    echo "[TEST]  FASTQC      = ${FASTQC}"
    echo "[TEST]  THREADS     = ${THREADS}"
fi

if ! [[ "${THREADS}" =~ ^[0-9]+$ ]] || [[ "${THREADS}" -lt 1 ]]; then
    echo "[ERROR] --threads must be a positive integer: ${THREADS}"
    exit 1
fi

if [[ ! -f "${PARAMS_YAML}" ]]; then
    echo "[ERROR] cellranger params YAML not found: ${PARAMS_YAML}"
    exit 1
fi

# ---- Output paths -----------------------------------------------------------
FASTQC_OUTDIR="${OUTDIR}/fastqc"
MULTIQC_OUTDIR="${OUTDIR}/multiqc"
MANIFEST_PATH="${OUTDIR}/fastq_manifest.tsv"
SUMMARY_PATH="${OUTDIR}/fastq_validation_summary.tsv"

# ---- Helpers ----------------------------------------------------------------
require_cmd() {
    local cmd="$1"
    local hint="$2"
    if ! command -v "${cmd}" &>/dev/null; then
        echo "[ERROR] Required command not found: ${cmd}"
        echo "        ${hint}"
        return 1
    fi
    return 0
}

# Infer read direction from common Illumina naming patterns so the summary can
# flag incomplete or mismatched paired-end inputs before Cell Ranger runs.
read_type_from_file() {
    local basename="$1"
    case "${basename}" in
        *_R1_*|*_R1.*|*_1.fastq.gz|*_1.fq.gz|*.R1.fastq.gz|*.R1.fq.gz)
            echo "R1"
            ;;
        *_R2_*|*_R2.*|*_2.fastq.gz|*_2.fq.gz|*.R2.fastq.gz|*.R2.fq.gz)
            echo "R2"
            ;;
        *)
            echo ""
            ;;
    esac
}

# Use Python for CSV parsing so we can validate column structure, tolerate
# UTF-8 BOMs, and emit one manifest row per declared sample/library pair.
build_manifest() {
    python3 - "$SAMPLES_CSV" <<'PY'
import csv
import sys

samples_csv = sys.argv[1]
required = ["sample_id", "gex_fastq_path", "gex_fastq_prefix"]
lib_specs = [
    ("GEX", "gex_fastq_path", "gex_fastq_prefix", True),
    ("VDJ-T", "vdj_t_fastq_path", "vdj_t_fastq_prefix", False),
    ("VDJ-T-GD", "vdj_t_gd_fastq_path", "vdj_t_gd_fastq_prefix", False),
    ("VDJ-B", "vdj_b_fastq_path", "vdj_b_fastq_prefix", False),
]

try:
    with open(samples_csv, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        missing = [name for name in required if name not in fields]
        if missing:
            print(
                "[ERROR] samples.csv missing required columns: "
                + ", ".join(missing),
                file=sys.stderr,
            )
            sys.exit(1)

        print("sample_id\tlibrary_type\tfastq_dir\tprefix")
        seen_ids = set()
        row_count = 0
        for line_no, row in enumerate(reader, start=2):
            row_count += 1
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id:
                print(
                    f"[ERROR] samples.csv line {line_no}: sample_id is empty",
                    file=sys.stderr,
                )
                sys.exit(1)
            if sample_id in seen_ids:
                print(
                    f"[ERROR] samples.csv line {line_no}: duplicate sample_id "
                    f"{sample_id}",
                    file=sys.stderr,
                )
                sys.exit(1)
            seen_ids.add(sample_id)

            for lib_type, path_col, prefix_col, required_pair in lib_specs:
                fastq_dir = (row.get(path_col) or "").strip()
                prefix = (row.get(prefix_col) or "").strip()
                if required_pair and (not fastq_dir or not prefix):
                    print(
                        f"[ERROR] samples.csv line {line_no}: {lib_type} needs "
                        f"both {path_col} and {prefix_col}",
                        file=sys.stderr,
                    )
                    sys.exit(1)
                if bool(fastq_dir) != bool(prefix):
                    print(
                        f"[ERROR] samples.csv line {line_no}: {lib_type} must "
                        f"set both {path_col} and {prefix_col} or leave both empty",
                        file=sys.stderr,
                    )
                    sys.exit(1)
                if fastq_dir and prefix:
                    print(
                        f"{sample_id}\t{lib_type}\t{fastq_dir}\t{prefix}"
                    )

        if row_count == 0:
            print("[ERROR] samples.csv has no sample rows", file=sys.stderr)
            sys.exit(1)
except FileNotFoundError:
    print(f"[ERROR] samples.csv not found: {samples_csv}", file=sys.stderr)
    sys.exit(1)
PY
}

# ---- Dry run ----------------------------------------------------------------
if [[ "${TEST_RUN}" == true ]]; then
    echo "[DRY]   Would validate FASTQ inputs from ${SAMPLES_CSV}"
    echo "[DRY]   Would read FastQC path from ${PARAMS_YAML}"
    echo "[DRY]   Would write reports to ${OUTDIR}"
    echo "[DRY]   Would run FastQC (${FASTQC}) with ${THREADS} thread(s)"
    echo "[DRY]   Would run MultiQC via: conda run -n multiqc multiqc"
    exit 0
fi

# ---- Dependency checks ------------------------------------------------------
ERRORS=0
WARNINGS=0

if ! require_cmd "python3" "Install python3 before running this script."; then
    ERRORS=$((ERRORS + 1))
fi
if [[ "${VALIDATE_ONLY}" == false ]]; then
    if [[ "${FASTQC}" == "${DEFAULT_FASTQC}" ]]; then
        echo "[WARN]  fastqc_path not set in ${PARAMS_YAML}; falling back to '${DEFAULT_FASTQC}'."
        WARNINGS=$((WARNINGS + 1))
    fi
    if [[ ! -x "${FASTQC}" ]] && ! command -v "${FASTQC}" &>/dev/null; then
        echo "[ERROR] FastQC not found: ${FASTQC}"
        echo "        Set cellranger.fastqc_path in ${PARAMS_YAML}, e.g.:"
        echo "        cellranger:"
        echo "          fastqc_path: /path/to/FastQC/fastqc"
        ERRORS=$((ERRORS + 1))
    fi
    if ! command -v "conda" &>/dev/null; then
        echo "[ERROR] conda not found - required to run MultiQC."
        echo "        Install Miniconda or Anaconda and ensure conda is on PATH."
        ERRORS=$((ERRORS + 1))
    elif ! conda run -n multiqc multiqc --version &>/dev/null 2>&1; then
        echo "[ERROR] conda environment 'multiqc' not found or multiqc not installed in it."
        echo "        Create with: conda create -n multiqc -c bioconda -c conda-forge multiqc"
        ERRORS=$((ERRORS + 1))
    fi
fi

[[ "${ERRORS}" -gt 0 ]] && exit 1

# ---- Initialize run artifacts ----------------------------------------------
mkdir -p "${OUTDIR}"
TMP_MANIFEST=$(mktemp)
trap 'rm -f "${TMP_MANIFEST}"' EXIT

if ! build_manifest > "${TMP_MANIFEST}"; then
    exit 1
fi

{
    printf "sample_id\tlibrary_type\tfastq_dir\tprefix\tfile_count\tr1_count\tr2_count\tstatus\n"
} > "${SUMMARY_PATH}"

{
    printf "sample_id\tlibrary_type\tfastq_dir\tprefix\tread_type\tfastq_file\n"
} > "${MANIFEST_PATH}"

# Validate each declared sample/library pair before any expensive FastQC work.
while IFS=$'\t' read -r sample_id library_type fastq_dir prefix; do
    if [[ "${sample_id}" == "sample_id" ]]; then
        continue
    fi

    status="OK"
    if [[ -z "${fastq_dir}" || -z "${prefix}" ]]; then
        echo "[ERROR] ${sample_id} ${library_type}: path/prefix unexpectedly empty"
        ERRORS=$((ERRORS + 1))
        printf "%s\t%s\t%s\t%s\t0\t0\t0\tERROR\n" \
            "${sample_id}" "${library_type}" "${fastq_dir}" "${prefix}" \
            >> "${SUMMARY_PATH}"
        continue
    fi

    if [[ "${fastq_dir}" == /path/to/* ]]; then
        echo "[ERROR] ${sample_id} ${library_type}: placeholder path still set: ${fastq_dir}"
        ERRORS=$((ERRORS + 1))
        printf "%s\t%s\t%s\t%s\t0\t0\t0\tERROR\n" \
            "${sample_id}" "${library_type}" "${fastq_dir}" "${prefix}" \
            >> "${SUMMARY_PATH}"
        continue
    fi

    if [[ ! -d "${fastq_dir}" ]]; then
        echo "[ERROR] ${sample_id} ${library_type}: FASTQ directory not found: ${fastq_dir}"
        ERRORS=$((ERRORS + 1))
        printf "%s\t%s\t%s\t%s\t0\t0\t0\tERROR\n" \
            "${sample_id}" "${library_type}" "${fastq_dir}" "${prefix}" \
            >> "${SUMMARY_PATH}"
        continue
    fi

    # Match both .fastq.gz and .fq.gz so the manifest follows common naming
    # conventions without requiring users to rename vendor output.
    shopt -s nullglob
    files=("${fastq_dir}/${prefix}"*.fastq.gz "${fastq_dir}/${prefix}"*.fq.gz)
    shopt -u nullglob

    if [[ "${#files[@]}" -eq 0 ]]; then
        echo "[ERROR] ${sample_id} ${library_type}: no FASTQ files matched prefix ${prefix} in ${fastq_dir}"
        ERRORS=$((ERRORS + 1))
        printf "%s\t%s\t%s\t%s\t0\t0\t0\tERROR\n" \
            "${sample_id}" "${library_type}" "${fastq_dir}" "${prefix}" \
            >> "${SUMMARY_PATH}"
        continue
    fi

    r1_count=0
    r2_count=0
    unknown_count=0

    for fastq_file in "${files[@]}"; do
        read_type=$(read_type_from_file "$(basename "${fastq_file}")")
        case "${read_type}" in
            R1) r1_count=$((r1_count + 1)) ;;
            R2) r2_count=$((r2_count + 1)) ;;
            *)
                unknown_count=$((unknown_count + 1))
                ;;
        esac

        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${sample_id}" "${library_type}" "${fastq_dir}" "${prefix}" \
            "${read_type:-UNKNOWN}" "${fastq_file}" \
            >> "${MANIFEST_PATH}"
    done

    if [[ "${r1_count}" -eq 0 || "${r2_count}" -eq 0 ]]; then
        echo "[ERROR] ${sample_id} ${library_type}: expected paired-end FASTQs with both R1 and R2"
        ERRORS=$((ERRORS + 1))
        status="ERROR"
    fi

    if [[ "${r1_count}" -ne "${r2_count}" ]]; then
        echo "[ERROR] ${sample_id} ${library_type}: R1/R2 file counts differ (${r1_count} vs ${r2_count})"
        ERRORS=$((ERRORS + 1))
        status="ERROR"
    fi

    if [[ "${unknown_count}" -gt 0 ]]; then
        echo "[WARN]  ${sample_id} ${library_type}: ${unknown_count} file(s) did not match R1/R2 naming heuristics"
        WARNINGS=$((WARNINGS + 1))
        if [[ "${status}" == "OK" ]]; then
            status="WARN"
        fi
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${sample_id}" "${library_type}" "${fastq_dir}" "${prefix}" \
        "${#files[@]}" "${r1_count}" "${r2_count}" "${status}" \
        >> "${SUMMARY_PATH}"

    # Only run FastQC after metadata and file matching checks pass for this sample/library pair.
    if [[ "${VALIDATE_ONLY}" == false && "${status}" != "ERROR" ]]; then
        mkdir -p "${FASTQC_OUTDIR}"
        echo "[RUN]   ${sample_id} ${library_type}: FastQC on ${#files[@]} file(s)"
        "${FASTQC}" --threads "${THREADS}" --outdir "${FASTQC_OUTDIR}" "${files[@]}"
    fi
done < "${TMP_MANIFEST}"

if [[ "${ERRORS}" -gt 0 ]]; then
    echo "[ERROR] FASTQ validation failed with ${ERRORS} error(s) and ${WARNINGS} warning(s)."
    echo "        Summary: ${SUMMARY_PATH}"
    exit 1
fi

if [[ "${VALIDATE_ONLY}" == true ]]; then
    echo "[OK]   FASTQ validation passed."
    echo "[OK]   Manifest: ${MANIFEST_PATH}"
    echo "[OK]   Summary:  ${SUMMARY_PATH}"
    exit 0
fi

# Aggregate per-file FastQC output into one report for the whole run.
mkdir -p "${MULTIQC_OUTDIR}"
echo "[RUN]   MultiQC summary"
conda run -n multiqc multiqc "${FASTQC_OUTDIR}" --outdir "${MULTIQC_OUTDIR}" --force

echo "[OK]   FASTQ QC complete."
echo "[OK]   Manifest: ${MANIFEST_PATH}"
echo "[OK]   Summary:  ${SUMMARY_PATH}"
echo "[OK]   FastQC:   ${FASTQC_OUTDIR}"
echo "[OK]   MultiQC:  ${MULTIQC_OUTDIR}"
