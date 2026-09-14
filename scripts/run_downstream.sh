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
targets_mode="false"
selected_step="${STEP:-}"
existing_run="${EXISTING_RUN:-}"
log_dir_explicit="false"
pyscenic_args=()

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/run_downstream.sh [options]

Options:
  --downstream-dir DIR   Downstream script directory, default: code/downstream
  --run-root DIR         Parent output directory, default: downstream_runs
  --run-dir DIR          Exact run output directory, default: <run-root>/<time>_pid<PID>
  --existing-run DIR     Reuse an existing run directory; accepts full path or name under --run-root
  --resume-run DIR       Alias for --existing-run
  --log-dir DIR          Log directory, default: <run-dir>/logs
  --rscript CMD          Rscript command, default: Rscript
  --email EMAIL          Final notification recipient, default: disabled
  --mail-cmd CMD         Mail command, default: msmtp
  --mail-from EMAIL      Sender address, default: MAIL_FROM or empty
  --step STEP            Run one step only; accepts 10, 10_MSH2, 10_MSH2.R, or a script path
  --targets              Build the compatibility DAG instead of running copied scripts
  --include-pyscenic     Also run 11_pySCENIC.sh and 12_pySCENIC_visualization.R
  --tf-list PATH         Passed to 11_pySCENIC.sh
  --ranking-db PATH      Passed to 11_pySCENIC.sh; may be repeated
  --ranking-dbs PATHS    Passed to 11_pySCENIC.sh
  --motif-annotations PATH
                         Passed to 11_pySCENIC.sh
  --project-dir DIR      Passed to 11_pySCENIC.sh; default for step 11 is <run-dir>
  --out-dir DIR          Passed to 11_pySCENIC.sh; default for step 11 is <run-dir>/pyscenic_10k
  --pyscenic-python PATH Passed to 11_pySCENIC.sh
  --pyscenic-command CMD Passed to 11_pySCENIC.sh
  --force-selection      Passed to 11_pySCENIC.sh to rebuild selected cells
  --force-export         Passed to 11_pySCENIC.sh to rebuild exported expression
  --pyscenic-arg ARG     Pass one extra argument to 11_pySCENIC.sh; may be repeated
  --continue-on-error    Run remaining scripts after failures, then send failure summary
  -h, --help             Show help

Examples:
  bash scripts/run_downstream.sh
  bash scripts/run_downstream.sh --step 10
  bash scripts/run_downstream.sh --step 10_MSH2.R
  bash scripts/run_downstream.sh --existing-run 20260825_213811_pid1045941 --step 10
  bash scripts/run_downstream.sh --step 11 --tf-list /path/to/allTFs_hg38.txt --ranking-db /path/to/ranking.feather --motif-annotations /path/to/motifs.tbl
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

make_abs_path() {
  local path="$1"
  case "$path" in
    /*) printf "%s\n" "$path" ;;
    .) printf "%s\n" "$script_dir" ;;
    ./*) printf "%s/%s\n" "$script_dir" "${path#./}" ;;
    *) printf "%s/%s\n" "$script_dir" "$path" ;;
  esac
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
  link_runtime_path "${script_dir}/code/downstream/lib" "${prepared_downstream_dir}/lib"
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
      --run-root)
        run_root="$2"
        run_dir="${2}/${run_id}"
        [[ "$log_dir_explicit" == "true" ]] || log_dir="${run_dir}/logs"
        shift 2
        ;;
      --run-dir)
        run_dir="$2"
        [[ "$log_dir_explicit" == "true" ]] || log_dir="${run_dir}/logs"
        shift 2
        ;;
      --existing-run|--resume-run)
        existing_run="$2"
        shift 2
        ;;
      --log-dir)
        log_dir="$2"
        log_dir_explicit="true"
        shift 2
        ;;
      --rscript) rscript_cmd="$2"; shift 2 ;;
      --email) notify_email="$2"; shift 2 ;;
      --mail-cmd) mail_cmd="$2"; shift 2 ;;
      --mail-from) mail_from="$2"; shift 2 ;;
      --step) selected_step="$2"; shift 2 ;;
      --targets) targets_mode="true"; shift ;;
      --include-pyscenic) include_pyscenic="true"; shift ;;
      --tf-list|--ranking-db|--ranking-dbs|--motif-annotations|--project-dir|--out-dir|--pyscenic-python|--pyscenic-command)
        [[ -n "${2:-}" ]] || die "Missing value for argument: $1"
        pyscenic_args+=("$1" "$2")
        shift 2
        ;;
      --force-selection|--force-export)
        pyscenic_args+=("$1")
        shift
        ;;
      --pyscenic-arg)
        [[ -n "${2:-}" ]] || die "Missing value for argument: $1"
        pyscenic_args+=("$2")
        shift 2
        ;;
      --continue-on-error) continue_on_error="true"; shift ;;
      -h|--help) usage; exit 0 ;;
      *) die "Unknown argument: $1" ;;
    esac
  done
}

resolve_existing_run_dir() {
  local requested="$1"

  [[ -n "$requested" ]] || return 1
  if [[ "$requested" = /* && -d "$requested" ]]; then
    make_abs_path "$requested"
    return 0
  fi
  if [[ -d "$requested" ]]; then
    make_abs_path "$requested"
    return 0
  fi
  if [[ -d "${run_root}/${requested}" ]]; then
    make_abs_path "${run_root}/${requested}"
    return 0
  fi

  return 1
}

has_pyscenic_arg() {
  local flag="$1"
  local arg
  for arg in "${pyscenic_args[@]}"; do
    [[ "$arg" == "$flag" ]] && return 0
  done
  return 1
}

build_pyscenic_run_args() {
  local args=()
  if ! has_pyscenic_arg "--project-dir"; then
    args+=("--project-dir" "$run_dir")
  fi
  if ! has_pyscenic_arg "--out-dir"; then
    args+=("--out-dir" "pyscenic_10k")
  fi
  args+=("${pyscenic_args[@]}")
  printf "%s\0" "${args[@]}"
}

build_pyscenic_visualization_args() {
  local args=()
  args+=("--project-dir" "$run_dir")
  args+=("--auc-loom" "pyscenic_10k/auc_mtx.loom")
  args+=("--regulons-csv" "pyscenic_10k/regulons.csv")
  args+=("--adj-tsv" "pyscenic_10k/adjacencies.tsv")
  printf "%s\0" "${args[@]}"
}

run_one_script() {
  local script="$1"
  local base log_file start_time end_time exit_code
  local args=()

  require_cmd "$rscript_cmd"
  base="$(basename "$script")"
  if [[ "$base" == "12_pySCENIC_visualization.R" ]]; then
    mapfile -d '' -t args < <(build_pyscenic_visualization_args)
  fi
  log_file="${log_dir}/${base}.log"
  start_time="$(date '+%F %T')"
  log "Running ${script}; log: ${log_file}"

  set +e
  {
    printf "[START] %s\n" "$start_time"
    printf "[CMD] %s %s" "$rscript_cmd" "$script"
    printf " %q" "${args[@]}"
    printf "\n"
    "$rscript_cmd" "$script" "${args[@]}"
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
  local args=()

  base="$(basename "$script")"
  if [[ "$base" == "11_pySCENIC.sh" ]]; then
    mapfile -d '' -t args < <(build_pyscenic_run_args)
  fi
  log_file="${log_dir}/${base}.log"
  start_time="$(date '+%F %T')"
  log "Running ${script}; log: ${log_file}"

  set +e
  {
    printf "[START] %s\n" "$start_time"
    printf "[CMD] bash %s" "$script"
    printf " %q" "${args[@]}"
    printf "\n"
    bash "$script" "${args[@]}"
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

run_script_by_type() {
  local script="$1"
  case "$script" in
    *.R) run_one_script "$script" ;;
    *.sh) run_pyscenic_script "$script" ;;
    *) die "Unsupported downstream script type: $script" ;;
  esac
}

resolve_selected_step() {
  local step="$1"
  local candidate candidate_base candidate_stem base stem number

  [[ -n "$step" ]] || return 1
  base="$(basename "$step")"
  stem="${base%.*}"

  for candidate in "${all_scripts[@]}"; do
    [[ -f "$candidate" ]] || continue
    candidate_base="$(basename "$candidate")"
    candidate_stem="${candidate_base%.*}"
    number="${candidate_base%%_*}"
    case "$step" in
      "$candidate"|"$candidate_base"|"$candidate_stem"|"$number")
        printf "%s\n" "$candidate"
        return 0
        ;;
    esac
    case "$base" in
      "$candidate_base"|"$candidate_stem"|"$number")
        printf "%s\n" "$candidate"
        return 0
        ;;
    esac
    case "$stem" in
      "$candidate_stem"|"$number")
        printf "%s\n" "$candidate"
        return 0
        ;;
    esac
  done

  return 1
}

parse_args "$@"
if [[ "$targets_mode" == "true" ]]; then
  [[ "$include_pyscenic" == "false" ]] || die "--targets does not yet manage pySCENIC; use the existing --include-pyscenic runner."
  [[ -z "$existing_run" ]] || die "--targets uses the configured project artifact paths and cannot be combined with --existing-run."
  target_args=()
  if [[ -n "$selected_step" ]]; then
    case "${selected_step%%_*}" in
      1) target_name="step_01_read_data" ;;
      2) target_name="step_02_data_clean" ;;
      3) target_name="step_03_cell_annotation" ;;
      4) target_name="step_04_vdj_gene" ;;
      5) target_name="step_05_cdr3_stat" ;;
      6) target_name="step_06_cdr3_between_samples" ;;
      7) target_name="step_07_cdr3_paired_sankey" ;;
      8) target_name="step_08_vd1_vs_vd2" ;;
      9) target_name="step_09_cdr3_paired" ;;
      10) target_name="step_10_msh2" ;;
      13) target_name="step_13_vd2_pseudotime" ;;
      14) target_name="step_14_vd1_vd2_pairwise" ;;
      15) target_name="step_15_cd80_cd86" ;;
      16) target_name="step_16_vd1_vd2_extra" ;;
      17) target_name="step_17_zol_pan_effector_vd2" ;;
      *) die "No targets compatibility mapping for downstream step: $selected_step" ;;
    esac
    target_args=(--target "$target_name")
  fi
  exec Rscript --vanilla "${script_dir}/code/pipeline/run.R" "${target_args[@]}"
fi
if [[ -n "$existing_run" ]]; then
  run_dir="$(resolve_existing_run_dir "$existing_run")" || die "Existing run directory not found: $existing_run"
  run_id="$(basename "$run_dir")"
  [[ "$log_dir_explicit" == "true" ]] || log_dir="${run_dir}/logs"
fi
run_dir="$(make_abs_path "$run_dir")"
if [[ "$log_dir_explicit" == "true" ]]; then
  log_dir="$(make_abs_path "$log_dir")"
else
  log_dir="${run_dir}/logs"
fi
prepared_downstream_dir="${run_dir}/downstream_scripts"
require_cmd perl
[[ -d "$downstream_dir" ]] || die "Downstream directory not found: $downstream_dir"
prepare_run_directory

summary_file="${run_dir}/run_summary.tsv"
message_file="${run_dir}/run_message.txt"
printf "script\tstatus\tstart_time\tend_time\tlog_file\n" > "$summary_file"

main_scripts=(
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
all_scripts=(
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
  "${prepared_downstream_dir}/11_pySCENIC.sh"
  "${prepared_downstream_dir}/12_pySCENIC_visualization.R"
  "${prepared_downstream_dir}/13_Vd2_pseudotime.R"
  "${prepared_downstream_dir}/14_Vd1Vd2_pairwise.R"
  "${prepared_downstream_dir}/15_CD80_CD86_expression.R"
  "${prepared_downstream_dir}/16_Vd1Vd2_extra_visualization.R"
  "${prepared_downstream_dir}/17_ZOL_PAN_effector_Vd2_comparison.R"
)

if [[ -n "$selected_step" ]]; then
  selected_script="$(resolve_selected_step "$selected_step")" || die "Unknown downstream step: $selected_step"
  scripts=("$selected_script")
  include_pyscenic="false"
else
  scripts=("${main_scripts[@]}")
fi

{
  printf "run_id\t%s\n" "$run_id"
  printf "run_dir\t%s\n" "$run_dir"
  printf "log_dir\t%s\n" "$log_dir"
  printf "start_time\t%s\n" "$(date '+%F %T')"
  printf "workdir\t%s\n" "$PWD"
  printf "downstream_dir\t%s\n" "$downstream_dir"
  printf "selected_step\t%s\n" "${selected_step:-all}"
  printf "include_pyscenic\t%s\n" "$include_pyscenic"
  printf "pyscenic_args\t"
  printf "%q " "${pyscenic_args[@]}"
  printf "\n"
  printf "continue_on_error\t%s\n" "$continue_on_error"
} > "${run_dir}/run_info.tsv"

failed=0
failed_script=""
failed_log=""

for script in "${scripts[@]}"; do
  [[ -f "$script" ]] || die "Missing downstream script: $script"
  if ! run_script_by_type "$script"; then
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
    if [[ -n "$selected_step" ]]; then
      printf "Selected downstream step completed successfully: %s\n" "$selected_step"
    else
      printf "All downstream scripts completed successfully.\n"
    fi
    printf "Run directory: %s\n" "$run_dir"
    printf "Logs: %s\n" "$log_dir"
    printf "Summary: %s\n" "$summary_file"
  } > "$message_file"
  send_email_notification "SUCCESS" "$message_file"
  if [[ -n "$selected_step" ]]; then
    log "Selected downstream step completed successfully: ${selected_step}. Summary: ${summary_file}"
  else
    log "All downstream scripts completed successfully. Summary: ${summary_file}"
  fi
else
  {
    if [[ -n "$selected_step" ]]; then
      printf "Selected downstream step failed: %s\n" "$selected_step"
    else
      printf "Downstream run failed.\n"
    fi
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
