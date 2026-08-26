#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$script_dir"

downstream_dir="${DOWNSTREAM_DIR:-code/downstream}"
run_root="${RUN_ROOT:-downstream_runs}"
run_id="${RUN_ID:-$(date '+%Y%m%d_%H%M%S')_pid$$}"
run_dir="${RUN_DIR:-${run_root}/${run_id}}"
log_dir="${LOG_DIR:-${run_dir}/logs}"
prepared_downstream_dir="${run_dir}/downstream_scripts"
rscript_cmd="${RSCRIPT_CMD:-Rscript}"
notify_email="${NOTIFY_EMAIL:-}"
mail_cmd="${MAIL_CMD:-msmtp}"
mail_from="${MAIL_FROM:-}"
include_pyscenic="${INCLUDE_PYSCENIC:-false}"
continue_on_error="${CONTINUE_ON_ERROR:-false}"

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/run_downstream.sh [options]

Options:
  --downstream-dir DIR   Downstream script directory, default: code/downstream
  --run-root DIR         Parent output directory, default: downstream_runs
  --run-dir DIR          Exact run output directory, default: <run-root>/<time>_pid<PID>
  --log-dir DIR          Log directory, default: <run-dir>/logs
  --rscript CMD          Rscript command, default: Rscript
  --email EMAIL          Final notification recipient, default: disabled
  --mail-cmd CMD         Mail command, default: msmtp
  --mail-from EMAIL      Sender address, default: MAIL_FROM or empty
  --include-pyscenic     Also run 11_pySCENIC.sh and 12_pySCENIC_visualization.R
  --continue-on-error    Run remaining scripts after failures, then send failure summary
  -h, --help             Show help

Examples:
  bash scripts/run_downstream.sh
  bash scripts/run_downstream.sh --include-pyscenic
  bash scripts/run_downstream.sh --run-root downstream_runs
  NOTIFY_EMAIL='user@example.com' MAIL_FROM='sender@example.com' bash scripts/run_downstream.sh
USAGE
}

log() {
  echo "[$(date '+%F %T')] $*"
}

die() {
  echo "[ERROR] $*" >&2
  exit 1
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

send_email_notification() {
  local status="$1"
  local detail_file="$2"
  local subject

  [[ -n "$notify_email" ]] || return 0
  if ! command -v "$mail_cmd" >/dev/null 2>&1; then
    echo "[WARN] Mail command not found, skipping email notification: ${mail_cmd}" >&2
    return 0
  fi

  subject="[run_downstream.sh] ${status}: $(basename "$PWD")"
  {
    printf "From: %s\n" "$mail_from"
    printf "To: %s\n" "$notify_email"
    printf "Subject: %s\n" "$subject"
    printf "\n"
    printf "Status: %s\n" "$status"
    printf "Time: %s\n" "$(date '+%F %T')"
    printf "Host: %s\n" "$(hostname)"
    printf "Workdir: %s\n" "$PWD"
    printf "Run directory: %s\n" "$run_dir"
    printf "Log directory: %s\n" "$log_dir"
    printf "\n"
    cat "$detail_file"
  } | "$mail_cmd" "$notify_email" || echo "[WARN] Failed to send email notification to ${notify_email}" >&2
}

link_runtime_path() {
  local source_path="$1"
  local link_path="$2"

  [[ -e "$link_path" || -L "$link_path" ]] && return 0
  ln -s "$source_path" "$link_path"
}

prepare_run_directory() {
  mkdir -p "$run_dir" "$log_dir" "$prepared_downstream_dir"
  link_runtime_path "${script_dir}/code" "${run_dir}/code"
  link_runtime_path "${script_dir}/code/config" "${run_dir}/config"
  link_runtime_path "${script_dir}/renv" "${run_dir}/renv"
  if [[ -f "${script_dir}/renv.lock" && ! -e "${run_dir}/renv.lock" ]]; then
    ln -s "${script_dir}/renv.lock" "${run_dir}/renv.lock"
  fi

  find "$downstream_dir" -maxdepth 1 \( -name "*.R" -o -name "*.sh" \) -type f -print | while IFS= read -r source_script; do
    local dest_script
    dest_script="${prepared_downstream_dir}/$(basename "$source_script")"
    cp "$source_script" "$dest_script"
    RUN_DIR_FOR_R="$run_dir" perl -0pi -e '
      my $d = $ENV{"RUN_DIR_FOR_R"};
      $d =~ s/\\/\\\\/g;
      $d =~ s/"/\\"/g;
      s#setwd\("/path/to/project"\)#"setwd(\"$d\")"#eg;
      s#project_dir <- "/path/to/project"#"project_dir <- \"$d\""#eg;
      s#project_dir = "/path/to/project"#"project_dir = \"$d\""#eg;
      s#PROJECT_DIR="\$\{PROJECT_DIR:-/path/to/project\}"#"PROJECT_DIR=\"\${PROJECT_DIR:-$d}\""#eg;
    ' "$dest_script"
    chmod --reference="$source_script" "$dest_script" 2>/dev/null || true
  done
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --downstream-dir) downstream_dir="$2"; shift 2 ;;
      --run-root) run_root="$2"; run_dir="${2}/${run_id}"; log_dir="${run_dir}/logs"; shift 2 ;;
      --run-dir) run_dir="$2"; log_dir="${run_dir}/logs"; shift 2 ;;
      --log-dir) log_dir="$2"; shift 2 ;;
      --rscript) rscript_cmd="$2"; shift 2 ;;
      --email) notify_email="$2"; shift 2 ;;
      --mail-cmd) mail_cmd="$2"; shift 2 ;;
      --mail-from) mail_from="$2"; shift 2 ;;
      --include-pyscenic) include_pyscenic="true"; shift ;;
      --continue-on-error) continue_on_error="true"; shift ;;
      -h|--help) usage; exit 0 ;;
      *) die "Unknown argument: $1" ;;
    esac
  done
}

run_one_script() {
  local script="$1"
  local base log_file start_time end_time exit_code

  base="$(basename "$script")"
  log_file="${log_dir}/${base}.log"
  start_time="$(date '+%F %T')"
  log "Running ${script}; log: ${log_file}"

  set +e
  {
    printf "[START] %s\n" "$start_time"
    printf "[CMD] %s %s\n" "$rscript_cmd" "$script"
    "$rscript_cmd" "$script"
  } > "$log_file" 2>&1
  exit_code="$?"
  set -e
  end_time="$(date '+%F %T')"

  if [[ "$exit_code" -eq 0 ]]; then
    printf "[END] %s\n[STATUS] SUCCESS\n" "$end_time" >> "$log_file"
    printf "%s\tSUCCESS\t%s\t%s\t%s\n" "$script" "$start_time" "$end_time" "$log_file" >> "$summary_file"
    log "Finished ${script}"
  else
    printf "[END] %s\n[STATUS] FAILED exit_code=%s\n" "$end_time" "$exit_code" >> "$log_file"
    printf "%s\tFAILED\t%s\t%s\t%s\n" "$script" "$start_time" "$end_time" "$log_file" >> "$summary_file"
    log "Failed ${script}; log: ${log_file}"
    return "$exit_code"
  fi
}

run_pyscenic_script() {
  local script="$1"
  local base log_file start_time end_time exit_code

  base="$(basename "$script")"
  log_file="${log_dir}/${base}.log"
  start_time="$(date '+%F %T')"
  log "Running ${script}; log: ${log_file}"

  set +e
  {
    printf "[START] %s\n" "$start_time"
    printf "[CMD] bash %s\n" "$script"
    bash "$script"
  } > "$log_file" 2>&1
  exit_code="$?"
  set -e
  end_time="$(date '+%F %T')"

  if [[ "$exit_code" -eq 0 ]]; then
    printf "[END] %s\n[STATUS] SUCCESS\n" "$end_time" >> "$log_file"
    printf "%s\tSUCCESS\t%s\t%s\t%s\n" "$script" "$start_time" "$end_time" "$log_file" >> "$summary_file"
    log "Finished ${script}"
  else
    printf "[END] %s\n[STATUS] FAILED exit_code=%s\n" "$end_time" "$exit_code" >> "$log_file"
    printf "%s\tFAILED\t%s\t%s\t%s\n" "$script" "$start_time" "$end_time" "$log_file" >> "$summary_file"
    log "Failed ${script}; log: ${log_file}"
    return "$exit_code"
  fi
}

parse_args "$@"
require_cmd "$rscript_cmd"
require_cmd perl
[[ -d "$downstream_dir" ]] || die "Downstream directory not found: $downstream_dir"
prepare_run_directory

summary_file="${run_dir}/run_summary.tsv"
message_file="${run_dir}/run_message.txt"
printf "script\tstatus\tstart_time\tend_time\tlog_file\n" > "$summary_file"
{
  printf "run_id\t%s\n" "$run_id"
  printf "run_dir\t%s\n" "$run_dir"
  printf "log_dir\t%s\n" "$log_dir"
  printf "start_time\t%s\n" "$(date '+%F %T')"
  printf "workdir\t%s\n" "$PWD"
  printf "downstream_dir\t%s\n" "$downstream_dir"
  printf "include_pyscenic\t%s\n" "$include_pyscenic"
  printf "continue_on_error\t%s\n" "$continue_on_error"
} > "${run_dir}/run_info.tsv"

scripts=(
  "${prepared_downstream_dir}/1_ReadData.R"
  "${prepared_downstream_dir}/2_DataClean.R"
  "${prepared_downstream_dir}/3_cellAnnotation.R"
  "${prepared_downstream_dir}/4_VDJgene.R"
  "${prepared_downstream_dir}/5_CDR3stat.R"
  "${prepared_downstream_dir}/6_CDR3betweenSamples.R"
  "${prepared_downstream_dir}/7_CDR3pairedSankey.R"
  "${prepared_downstream_dir}/8_Vd1vs2.R"
  "${prepared_downstream_dir}/9_CDR3paired.R"
  "${prepared_downstream_dir}/10_MSH2.R"
  "${prepared_downstream_dir}/13_Vd2_pseudotime.R"
  "${prepared_downstream_dir}/14_Vd1Vd2_pairwise.R"
  "${prepared_downstream_dir}/15_CD80_CD86_expression.R"
  "${prepared_downstream_dir}/16_Vd1Vd2_extra_visualization.R"
  "${prepared_downstream_dir}/17_ZOL_PAN_effector_Vd2_comparison.R"
)

failed=0
failed_script=""
failed_log=""

for script in "${scripts[@]}"; do
  [[ -f "$script" ]] || die "Missing downstream script: $script"
  if ! run_one_script "$script"; then
    failed=1
    failed_script="$script"
    failed_log="${log_dir}/$(basename "$script").log"
    [[ "$continue_on_error" == "true" ]] || break
  fi
done

if [[ "$failed" -eq 0 && "$include_pyscenic" == "true" ]]; then
  pyscenic_script="${prepared_downstream_dir}/11_pySCENIC.sh"
  [[ -f "$pyscenic_script" ]] || die "Missing pySCENIC script: $pyscenic_script"
  if ! run_pyscenic_script "$pyscenic_script"; then
    failed=1
    failed_script="$pyscenic_script"
    failed_log="${log_dir}/$(basename "$pyscenic_script").log"
  fi
fi

if [[ "$failed" -eq 0 && "$include_pyscenic" == "true" ]]; then
  pyscenic_vis_script="${prepared_downstream_dir}/12_pySCENIC_visualization.R"
  [[ -f "$pyscenic_vis_script" ]] || die "Missing pySCENIC visualization script: $pyscenic_vis_script"
  if ! run_one_script "$pyscenic_vis_script"; then
    failed=1
    failed_script="$pyscenic_vis_script"
    failed_log="${log_dir}/$(basename "$pyscenic_vis_script").log"
  fi
fi

if [[ "$failed" -eq 0 ]]; then
  {
    printf "All downstream scripts completed successfully.\n"
    printf "Run directory: %s\n" "$run_dir"
    printf "Logs: %s\n" "$log_dir"
    printf "Summary: %s\n" "$summary_file"
  } > "$message_file"
  send_email_notification "SUCCESS" "$message_file"
  log "All downstream scripts completed successfully. Summary: ${summary_file}"
else
  {
    printf "Downstream run failed.\n"
    printf "Run directory: %s\n" "$run_dir"
    printf "Logs: %s\n" "$log_dir"
    printf "Failed script: %s\n" "$failed_script"
    printf "Failed log: %s\n" "$failed_log"
    printf "Summary: %s\n" "$summary_file"
    printf "\nLast 60 log lines:\n"
    tail -n 60 "$failed_log"
  } > "$message_file"
  send_email_notification "FAILED" "$message_file"
  exit 1
fi
