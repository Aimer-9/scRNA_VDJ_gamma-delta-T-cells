#!/bin/bash
# Run cellranger multi sequentially for each generated sample config.

set -uo pipefail

PARAMS_YAML="config/params.yaml"
CONFIG_DIR=""
OUTDIR_ROOT_ARG=""
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --params) PARAMS_YAML="$2"; shift 2 ;;
        --configs) CONFIG_DIR="$2"; shift 2 ;;
        --outdir) OUTDIR_ROOT_ARG="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help)
            cat <<'USAGE'
Usage:
  bash cellranger/src/03_cellranger_multi.sh [OPTIONS]

Options:
  --params PATH    Path to params.yaml [default: config/params.yaml]
  --configs DIR   Directory containing *_multi_config.csv files
  --outdir DIR    Pipeline output root
  --dry-run       Print planned commands without running Cell Ranger
  -h, --help      Show help
USAGE
            exit 0
            ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

yaml_get() {
    python3 - "$PARAMS_YAML" "$1" "$2" <<'PY'
import sys

params, dotted_key, default = sys.argv[1:]
try:
    import yaml
    with open(params, encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    value = data
    for key in dotted_key.split("."):
        if not isinstance(value, dict) or key not in value or value[key] is None:
            print(default)
            sys.exit(0)
        value = value[key]
    print(value)
except Exception:
    print(default)
PY
}

format_elapsed() {
    local seconds="$1"
    local hours=$((seconds / 3600))
    local minutes=$(((seconds % 3600) / 60))
    local secs=$((seconds % 60))
    if [[ "${hours}" -gt 0 ]]; then
        printf "%dh%02dm%02ds" "${hours}" "${minutes}" "${secs}"
    elif [[ "${minutes}" -gt 0 ]]; then
        printf "%dm%02ds" "${minutes}" "${secs}"
    else
        printf "%ds" "${secs}"
    fi
}

CELLRANGER="$(yaml_get cellranger.cellranger_path cellranger)"
LOCALCORES="$(yaml_get cellranger.localcores 16)"
LOCALMEM="$(yaml_get cellranger.localmem 64)"

if [[ -n "${OUTDIR_ROOT_ARG}" ]]; then
    OUTDIR_ROOT="${OUTDIR_ROOT_ARG}"
else
    OUTDIR_ROOT="$(yaml_get outdir output)"
fi

if [[ -z "${CONFIG_DIR}" ]]; then
    CONFIG_DIR="${OUTDIR_ROOT}/cellranger_configs"
fi

OUTPUT_DIR="${OUTDIR_ROOT}/cellranger_output"
LOG_DIR="${OUTPUT_DIR}/logs"
LOCK_FILE="${OUTPUT_DIR}/cellranger_done.lock"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"
OUTPUT_DIR="$(cd "${OUTPUT_DIR}" && pwd)"
LOG_DIR="${OUTPUT_DIR}/logs"
LOCK_FILE="${OUTPUT_DIR}/cellranger_done.lock"

if [[ ! -f "${LOCK_FILE}" ]]; then
    printf "sample_id\tcompleted_at\tstatus\tconfig_csv\n" > "${LOCK_FILE}"
fi

if [[ "${DRY_RUN}" == false ]]; then
    if ! command -v "${CELLRANGER}" >/dev/null 2>&1 && [[ ! -x "${CELLRANGER}" ]]; then
        echo "[ERROR] cellranger not found or not executable: ${CELLRANGER}" >&2
        echo "        Set cellranger.cellranger_path in ${PARAMS_YAML}" >&2
        exit 1
    fi
    if [[ ! -d "${CONFIG_DIR}" ]]; then
        echo "[ERROR] Config directory not found: ${CONFIG_DIR}" >&2
        exit 1
    fi
fi

shopt -s nullglob
configs=("${CONFIG_DIR}"/*_multi_config.csv)
shopt -u nullglob

if [[ "${#configs[@]}" -eq 0 ]]; then
    echo "[INFO]  No config CSVs found in ${CONFIG_DIR} — nothing to run."
    exit 0
fi

failures=()

for config in "${configs[@]}"; do
    sample_start_epoch="$(date +%s)"
    sample_id="$(basename "${config}" _multi_config.csv)"
    config_abs="$(realpath "${config}")"
    sample_output_dir="${OUTPUT_DIR}/${sample_id}"
    sample_log="${LOG_DIR}/cellranger_${sample_id}.log"

    if [[ "${DRY_RUN}" == false ]]; then
        if awk -F '\t' -v sid="${sample_id}" 'NR > 1 && $1 == sid && $3 == "DONE" { found = 1 } END { exit(found ? 0 : 1) }' "${LOCK_FILE}"; then
            elapsed="$(format_elapsed "$(( $(date +%s) - sample_start_epoch ))")"
            echo "[SKIP]  ${sample_id} — recorded complete (${elapsed})"
            {
                echo "[$(date -Is)] [SKIP] ${sample_id}"
                echo "Reason: sample is recorded as DONE in ${LOCK_FILE}."
                echo "Output: ${sample_output_dir}"
                echo "Elapsed: ${elapsed}"
            } >> "${sample_log}"
            continue
        fi

        if [[ -f "${sample_output_dir}/_SUCCESS" ]]; then
            elapsed="$(format_elapsed "$(( $(date +%s) - sample_start_epoch ))")"
            echo "[SKIP]  ${sample_id} — already complete (${elapsed})"
            {
                echo "[$(date -Is)] [SKIP] ${sample_id}"
                echo "Reason: existing Cell Ranger success marker found before launch."
                echo "Output: ${sample_output_dir}"
                echo "Elapsed: ${elapsed}"
            } >> "${sample_log}"
            printf "%s\t%s\tDONE\t%s\n" "${sample_id}" "$(date -Is)" "${config_abs}" >> "${LOCK_FILE}"
            continue
        fi

        if [[ -d "${sample_output_dir}/outs" ]]; then
            elapsed="$(format_elapsed "$(( $(date +%s) - sample_start_epoch ))")"
            echo "[SKIP]  ${sample_id} — output already exists (${elapsed})"
            {
                echo "[$(date -Is)] [SKIP] ${sample_id}"
                echo "Reason: existing Cell Ranger output directory found before launch."
                echo "Output: ${sample_output_dir}/outs"
                echo "Remove or rename ${sample_output_dir} to rerun this sample."
                echo "Elapsed: ${elapsed}"
            } >> "${sample_log}"
            continue
        fi
    fi

    if [[ "${DRY_RUN}" == true ]]; then
        echo "[DRY]   ${sample_id} — would run cellranger multi"
        echo "        log: ${sample_log}"
        continue
    fi

    echo "[RUN]   ${sample_id} — starting cellranger multi"

    (cd "${OUTPUT_DIR}" && "${CELLRANGER}" multi \
        --id="${sample_id}" \
        --csv="${config_abs}" \
        --localcores="${LOCALCORES}" \
        --localmem="${LOCALMEM}") > "${sample_log}" 2>&1
    status="$?"
    elapsed="$(format_elapsed "$(( $(date +%s) - sample_start_epoch ))")"

    if [[ "${status}" == "0" ]]; then
        echo "[DONE]  ${sample_id} — ${elapsed}"
        printf "%s\t%s\tDONE\t%s\n" "${sample_id}" "$(date -Is)" "${config_abs}" >> "${LOCK_FILE}"
    else
        echo "[FAIL]  ${sample_id} — status ${status}, ${elapsed}"
        failures+=("${sample_id}:${status}")
    fi
done

if [[ "${#failures[@]}" -gt 0 ]]; then
    echo "[ERROR] Cell Ranger failed for ${#failures[@]} sample(s): ${failures[*]}" >&2
    echo "        Logs are in: ${LOG_DIR}" >&2
    exit 1
fi

echo "=== All samples processed ==="
