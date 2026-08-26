#!/usr/bin/env bash
set -euo pipefail

# Integrate lane FASTQs by sample-library, then package and upload.
#
# Expected input:
#   sample-5LIB_S7_L004_R1_001.fastq.gz
#   sample-5LIB_S7_L004_R2_001.fastq.gz
#   sample-5LIB_S7_L004_I1_001.fastq.gz
#   sample-5LIB_S7_L004_I2_001.fastq.gz
#
# Integrated output:
#   ./integrated_fastq/sample-5LIB_R1_001.fastq.gz
#   ./integrated_fastq/sample-5LIB_R2_001.fastq.gz
#   ./integrated_fastq/sample-5LIB_I1_001.fastq.gz
#   ./integrated_fastq/sample-5LIB_I2_001.fastq.gz
#
# Package output:
#   ./upload_tar/sample-5LIB.tar.gz

config_file="${FTP_CONFIG:-$(dirname "${BASH_SOURCE[0]}")/env/ftp.env}"
for arg_index in "$@"; do
  if [[ "${previous_arg:-}" == "--config" ]]; then
    config_file="$arg_index"
    break
  fi
  previous_arg="$arg_index"
done
unset previous_arg || true

if [[ -f "$config_file" ]]; then
  # shellcheck source=/dev/null
  source "$config_file"
fi

rawdata_dir="${RAWDATA_DIR:-rawdata}"
result_dir="${RESULT_DIR:-.}"
integrated_dir="${INTEGRATED_DIR:-${result_dir}/integrated_fastq}"
upload_dir="${UPLOAD_DIR:-${result_dir}/upload_tar}"
jobs="${JOBS:-20}"

# FTP credentials are loaded from ftp.env by default. They can still be
# overridden by environment variables or CLI flags.
ftp_host="${FTP_HOST:-}"
ftp_user="${FTP_USER:-}"
ftp_pass="${FTP_PASS:-}"
remote_dir="${REMOTE_DIR:-GSA}"

# Email is optional at runtime. Set NOTIFY_EMAIL='' to suppress notifications,
# or override the recipient/sender with --email/--mail-from.
notify_email="${NOTIFY_EMAIL-}"
mail_cmd="${MAIL_CMD:-msmtp}"
mail_from="${MAIL_FROM:-${ftp_user}}"

tar_pe="${TAR_PE:-false}"
manifest="${upload_dir}/manifest.tsv"
sample_list="${result_dir}/sample_library.list"
pid_file="${result_dir}/ftp.sh.pid"

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/ftp.sh [integrate|tar|ftp|check|test-email|all|stop] [options]

Commands:
  integrate  Concatenate all lanes by sample-library and read type in parallel
  tar        Create md5 files and tar packages from integrated FASTQs
  ftp        Upload tar directory with lftp
  check      Check package structure
  test-email Send a test email notification
  all        integrate, tar, check, then ftp
  stop       Stop a running ftp.sh job and child processes

Options:
  --rawdata DIR       Input FASTQ directory, default: rawdata
  --result DIR        Working/output root, default: .
  --integrated DIR    Integrated FASTQ directory, default: <result>/integrated_fastq
  --upload DIR        Tar/upload directory, default: <result>/upload_tar
  --jobs INT          Parallel jobs, default: 20
  --tar-pe            Also tar R1/R2-only paired-end samples
  --config FILE       Source FTP/mail config file, default: scripts/env/ftp.env
  --ftp-host URL      FTP URL, default: FTP_HOST from config/env
  --ftp-user USER     FTP username
  --ftp-pass PASS     FTP password
  --remote-dir DIR    Remote FTP directory
  --email EMAIL       Send step result notifications
  --mail-cmd CMD      Mail command, default: msmtp
  --mail-from EMAIL   Sender address in email headers
  -h, --help          Show help

Examples:
  bash scripts/ftp.sh integrate --rawdata rawdata --result . --jobs 20
  bash scripts/ftp.sh tar --rawdata rawdata --result . --jobs 20
  bash scripts/ftp.sh check
  bash scripts/ftp.sh ftp
  bash scripts/ftp.sh test-email
  bash scripts/ftp.sh stop
  FTP_CONFIG=scripts/env/ftp.env REMOTE_DIR='GSA/subHRAxxxxxx' bash scripts/ftp.sh all
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

cleanup_pid_file() {
  if [[ -f "$pid_file" ]] && [[ "$(cat "$pid_file")" == "$$" ]]; then
    rm -f "$pid_file"
  fi
}

register_running_pid() {
  # Long-running commands write a PID file so accidental duplicate runs do not
  # write to the same integrated/upload directories concurrently.
  mkdir -p "$result_dir"
  if [[ -f "$pid_file" ]]; then
    local old_pid
    old_pid="$(cat "$pid_file")"
    if [[ -n "$old_pid" ]] && kill -0 "$old_pid" >/dev/null 2>&1; then
      die "Another ftp.sh run is active with PID ${old_pid}. Stop it with: bash scripts/ftp.sh stop"
    fi
  fi
  echo "$$" > "$pid_file"
  trap cleanup_pid_file EXIT
  trap 'echo "[STOP] Interrupted; stopping current run." >&2; exit 130' INT TERM
}

list_descendant_pids() {
  # GNU parallel and lftp can leave child processes behind. Walk the process
  # tree so "bash scripts/ftp.sh stop" can terminate the whole active run.
  local parent="$1"
  local child
  pgrep -P "$parent" 2>/dev/null | while IFS= read -r child; do
    list_descendant_pids "$child"
    echo "$child"
  done
}

stop_running_job() {
  require_cmd pgrep

  if [[ ! -f "$pid_file" ]]; then
    echo "[INFO] No running ftp.sh PID file found: ${pid_file}"
    return 0
  fi

  local pid descendants targets still_running
  pid="$(cat "$pid_file")"
  if [[ -z "$pid" ]] || ! kill -0 "$pid" >/dev/null 2>&1; then
    echo "[INFO] Stale PID file removed: ${pid_file}"
    rm -f "$pid_file"
    return 0
  fi

  descendants="$(list_descendant_pids "$pid" | sort -rn || true)"
  targets="$(printf "%s\n%s\n" "$descendants" "$pid" | awk 'NF')"

  echo "[STOP] Sending TERM to ftp.sh PID ${pid} and child processes."
  printf "%s\n" "$targets" | xargs -r kill -TERM
  sleep 3

  still_running="$(printf "%s\n" "$targets" | while IFS= read -r target; do
    if kill -0 "$target" >/dev/null 2>&1; then
      echo "$target"
    fi
  done)"
  if [[ -n "$still_running" ]]; then
    echo "[STOP] Sending KILL to remaining child processes."
    printf "%s\n" "$still_running" | xargs -r kill -KILL
  fi

  rm -f "$pid_file"
  echo "[DONE] Stopped ftp.sh run."
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --rawdata) rawdata_dir="$2"; shift 2 ;;
      --result) result_dir="$2"; integrated_dir="${2}/integrated_fastq"; upload_dir="${2}/upload_tar"; sample_list="${2}/sample_library.list"; manifest="${2}/upload_tar/manifest.tsv"; pid_file="${2}/ftp.sh.pid"; shift 2 ;;
      --integrated) integrated_dir="$2"; shift 2 ;;
      --upload) upload_dir="$2"; manifest="${2}/manifest.tsv"; shift 2 ;;
      --jobs) jobs="$2"; shift 2 ;;
      --tar-pe) tar_pe="true"; shift ;;
      --config) shift 2 ;;
      --ftp-host) ftp_host="$2"; shift 2 ;;
      --ftp-user) ftp_user="$2"; shift 2 ;;
      --ftp-pass) ftp_pass="$2"; shift 2 ;;
      --remote-dir) remote_dir="$2"; shift 2 ;;
      --email) notify_email="$2"; shift 2 ;;
      --mail-cmd) mail_cmd="$2"; shift 2 ;;
      --mail-from) mail_from="$2"; shift 2 ;;
      -h|--help) usage; exit 0 ;;
      *) die "Unknown argument: $1" ;;
    esac
  done
}

validate_jobs() {
  [[ "$jobs" =~ ^[0-9]+$ ]] && [[ "$jobs" -ge 1 ]] || die "--jobs must be a positive integer: $jobs"
}

sample_library_from_fastq() {
  # Convert 10x lane FASTQ names to sample-library IDs:
  #   EG-4-TCR_S101_L002_R1_001.fastq.gz -> EG-4-TCR
  sed -E 's/_S[0-9]+_L[0-9]+_[IR][0-9]+_001\.fastq\.gz$//'
}

read_type_from_fastq() {
  # Extract I1/I2/R1/R2 from a 10x FASTQ filename.
  sed -E 's/^.*_([IR][0-9]+)_001\.fastq\.gz$/\1/'
}
export -f read_type_from_fastq

make_sample_list() {
  # The tar step packages per sample-library, so build a stable unique sample
  # list from raw FASTQ symlinks/files before any integration starts.
  mkdir -p "$result_dir"
  find "$rawdata_dir" -maxdepth 1 \( -type f -o -type l \) -name "*.fastq.gz" -printf "%f\n" \
    | sample_library_from_fastq \
    | sort -u > "$sample_list"
  [[ -s "$sample_list" ]] || die "No FASTQ files found under: $rawdata_dir"
}

make_expected_integrated_list() {
  # For every raw lane FASTQ, predict the lane-collapsed integrated filename.
  # Multiple lanes of the same sample/read collapse to one expected output.
  find "$rawdata_dir" -maxdepth 1 \( -type f -o -type l \) -name "*.fastq.gz" -printf "%f\n" \
    | perl -ne 'chomp; if (/^(.+)_S\d+_L\d+_([IR]\d+)_001\.fastq\.gz$/) { print $1 . "_" . $2 . "_001.fastq.gz\n" } else { warn "[WARN] FASTQ name did not match 10x pattern: $_\n" }' \
    | sort -u
}

remove_incomplete_integrated_tmp_files() {
  # Integration writes to *.tmp.<pid> first. If a run is interrupted, those temp
  # files are discarded on the next resume so they cannot be mistaken for final
  # integrated FASTQs.
  local tmp_count

  [[ -d "$integrated_dir" ]] || return 0
  tmp_count="$(find "$integrated_dir" -maxdepth 1 -type f -name '*.fastq.gz.tmp.*' | wc -l)"
  if [[ "$tmp_count" -gt 0 ]]; then
    log "Removing ${tmp_count} incomplete temporary integrated FASTQ files before resume."
    find "$integrated_dir" -maxdepth 1 -type f -name '*.fastq.gz.tmp.*' -delete
  fi
}

check_integrated_expected_fastq() {
  # Worker used by GNU parallel. It emits one tab-delimited status line, which
  # lets the parent process summarize expected output presence quickly.
  local expected_fastq="$1"
  local integrated_dir="$2"
  local output_path="${integrated_dir}/${expected_fastq}"

  if [[ ! -s "$output_path" ]]; then
    printf "MISSING\t%s\n" "$output_path"
  else
    printf "COMPLETE\t%s\n" "$output_path"
  fi
}
export -f check_integrated_expected_fastq

audit_integrated_outputs() {
  # Compare expected integrated FASTQs against files on disk. In "resume" mode
  # missing files are only summarized; in "check" mode they fail the command.
  local mode="${1:-report}"
  local expected_list result_list result_status result_path complete missing status

  expected_list="$(mktemp)"
  result_list="$(mktemp)"
  make_expected_integrated_list > "$expected_list"

  if [[ -s "$expected_list" ]]; then
    parallel --halt never -j "$jobs" check_integrated_expected_fastq {} "$integrated_dir" :::: "$expected_list" > "$result_list"
  fi

  complete="$(awk -F '\t' '$1 == "COMPLETE" {n++} END {print n + 0}' "$result_list")"
  missing="$(awk -F '\t' '$1 == "MISSING" {n++} END {print n + 0}' "$result_list")"
  status=0

  while IFS=$'\t' read -r result_status result_path; do
    [[ -n "${result_status:-}" ]] || continue
    case "$result_status" in
      MISSING)
        [[ "$mode" == "check" ]] && echo "[FAIL] Missing integrated FASTQ: ${result_path}" >&2
        status=1
        ;;
    esac
  done < "$result_list"

  rm -f "$expected_list" "$result_list"
  echo "[INFO] Integrated audit: complete=${complete}, missing=${missing}"

  if [[ "$mode" == "check" ]]; then
    return "$status"
  fi
  return 0
}

integrate_one_read_type() {
  # Concatenate all lanes for one sample-library/read type. rawdata entries may
  # be symlinks; cat follows symlinks and writes a real integrated FASTQ.
  local sample_library="$1"
  local read_type="$2"
  local rawdata_dir="$3"
  local integrated_dir="$4"
  local out_file="${integrated_dir}/${sample_library}_${read_type}_001.fastq.gz"
  local tmp_file file_list

  file_list="$(find "$rawdata_dir" -maxdepth 1 \( -type f -o -type l \) \
    -name "${sample_library}_S*_L*_${read_type}_001.fastq.gz" -printf "%f\n" | sort)"

  [[ -n "$file_list" ]] || return 0

  mkdir -p "$integrated_dir"
  if [[ -s "$out_file" ]]; then
    # Resume behavior: existing non-empty final outputs are not rebuilt.
    echo "[SKIP] Complete integrated FASTQ exists: ${out_file}"
    return 0
  fi

  # Atomic-ish output: write temp, then rename to final. Interrupted runs leave
  # only *.tmp.<pid>, which are cleaned before resume.
  tmp_file="${out_file}.tmp.$$"
  rm -f "$tmp_file"
  while IFS= read -r fastq; do
    cat "${rawdata_dir}/${fastq}"
  done <<< "$file_list" > "$tmp_file"

  if [[ -s "$tmp_file" ]]; then
    mv -f "$tmp_file" "$out_file"
  else
    rm -f "$tmp_file"
    echo "[ERROR] Integrated FASTQ is empty: ${out_file}" >&2
    return 1
  fi

  echo "[OK] ${out_file}"
}
export -f integrate_one_read_type

integrate_one_sample() {
  # Detect which read types actually exist for a sample. Some project samples
  # have I1/I2/R1/R2, while PE-only samples have only R1/R2.
  local sample_library="$1"
  local rawdata_dir="$2"
  local integrated_dir="$3"
  local jobs="$4"
  local read_types

  read_types="$(find "$rawdata_dir" -maxdepth 1 \( -type f -o -type l \) \
    -name "${sample_library}_S*_L*_[IR]*_001.fastq.gz" -printf "%f\n" \
    | read_type_from_fastq | sort -u)"

  [[ -n "$read_types" ]] || {
    echo "[WARN] No read types found for ${sample_library}" >&2
    return 0
  }

  printf "%s\n" "$read_types" \
    | parallel --halt soon,fail=1 -j "$jobs" integrate_one_read_type "$sample_library" {} "$rawdata_dir" "$integrated_dir"
}
export -f integrate_one_sample

integrate_all() {
  require_cmd find
  require_cmd parallel
  require_cmd cat
  require_cmd perl
  require_cmd awk
  validate_jobs
  [[ -d "$rawdata_dir" ]] || die "Input FASTQ directory not found: $rawdata_dir"

  make_sample_list
  mkdir -p "$integrated_dir"
  remove_incomplete_integrated_tmp_files
  # Report current state before resuming, then rebuild only missing
  # sample/read outputs via integrate_one_read_type.
  audit_integrated_outputs "resume"

  log "Integrating lanes by sample-library and read type in parallel."
  parallel --halt soon,fail=1 -j "$jobs" integrate_one_sample {} "$rawdata_dir" "$integrated_dir" "$jobs" :::: "$sample_list"
  audit_integrated_outputs "check" || return 1
  log "Integrated FASTQs written to: $integrated_dir"
}

tar_one_sample() {
  # Package one sample-library. Four-read 10x samples become sample.tar.gz with
  # I1/I2/R1/R2 plus sample_md5.txt. PE-only samples are copied as plain FASTQs
  # unless --tar-pe is used.
  local sample_library="$1"
  local integrated_dir="$2"
  local upload_dir="$3"
  local tar_pe="$4"
  local manifest_tmp="$5"
  local tmpdir tar_file read_count fastq_count

  tmpdir="$(mktemp -d)"
  find "$integrated_dir" -maxdepth 1 -type f -name "${sample_library}_[IR]*_001.fastq.gz" -printf "%f\n" \
    | sort > "${tmpdir}/fastq.list"

  fastq_count="$(wc -l < "${tmpdir}/fastq.list")"
  if [[ "$fastq_count" -eq 0 ]]; then
    echo "[WARN] No integrated FASTQs found for ${sample_library}" >&2
    rm -rf "$tmpdir"
    return 0
  fi

  read_count="$(sed -E 's/^.*_([IR][0-9]+)_001\.fastq\.gz$/\1/' "${tmpdir}/fastq.list" | sort -u | wc -l)"
  if [[ "$read_count" -eq 2 && "$tar_pe" != "true" ]]; then
    # Reviewer note: plain paired-end R1/R2 data may not need compression.
    mkdir -p "${upload_dir}/plain_pe"
    while IFS= read -r fastq; do
      cp -f "${integrated_dir}/${fastq}" "${upload_dir}/plain_pe/${fastq}"
    done < "${tmpdir}/fastq.list"
    echo "[PE] ${sample_library}: copied R1/R2 to ${upload_dir}/plain_pe; not tarred."
    rm -rf "$tmpdir"
    return 0
  fi

  while IFS= read -r fastq; do
    # MD5 paths inside the md5 text are written as ./fastq_name, matching GSA
    # submission expectations after extracting the tar.
    md5sum "${integrated_dir}/${fastq}" | awk -v name="./${fastq}" '{print $1 "  " name}'
  done < "${tmpdir}/fastq.list" > "${tmpdir}/${sample_library}_md5.txt"

  tar_file="${upload_dir}/${sample_library}.tar.gz"
  mkdir -p "$upload_dir"
  tar -czf "$tar_file" \
    -C "$integrated_dir" $(tr '\n' ' ' < "${tmpdir}/fastq.list") \
    -C "$tmpdir" "${sample_library}_md5.txt"

  printf "%s\t%s\t%s\t%s\n" "$sample_library" "$(basename "$tar_file")" "$fastq_count" "$(md5sum "$tar_file" | awk '{print $1}')" >> "$manifest_tmp"
  echo "[TAR] ${tar_file}"
  rm -rf "$tmpdir"
}
export -f tar_one_sample

tar_all() {
  require_cmd find
  require_cmd parallel
  require_cmd md5sum
  require_cmd tar
  require_cmd perl
  require_cmd awk
  validate_jobs

  [[ -d "$rawdata_dir" ]] || die "Input FASTQ directory not found: $rawdata_dir"
  [[ -d "$integrated_dir" ]] || die "Integrated FASTQ directory not found: $integrated_dir"
  remove_incomplete_integrated_tmp_files
  # Do not package partial results. This catches stopped integrations before
  # tar files are created.
  audit_integrated_outputs "check" || return 1
  make_sample_list
  mkdir -p "$upload_dir"
  printf "sample_library\ttar_file\tfastq_count\ttar_md5\n" > "${manifest}.tmp"

  log "Creating MD5 files and tar packages in parallel."
  parallel --halt soon,fail=1 -j "$jobs" tar_one_sample {} "$integrated_dir" "$upload_dir" "$tar_pe" "${manifest}.tmp" :::: "$sample_list"
  mv "${manifest}.tmp" "$manifest"
  log "Package manifest written to: $manifest"
}

check_all() {
  require_cmd tar
  require_cmd perl
  require_cmd parallel
  require_cmd awk
  log "Checking outputs."

  [[ -d "$rawdata_dir" ]] || die "Input FASTQ directory not found: $rawdata_dir"
  [[ -d "$integrated_dir" ]] || die "Missing integrated directory: $integrated_dir"
  [[ -d "$upload_dir" ]] || die "Missing upload directory: $upload_dir"

  echo "[INFO] Integrated FASTQ count: $(find "$integrated_dir" -maxdepth 1 -type f -name '*.fastq.gz' | wc -l)"
  local tmp_count
  tmp_count="$(find "$integrated_dir" -maxdepth 1 -type f -name '*.fastq.gz.tmp.*' | wc -l)"
  if [[ "$tmp_count" -gt 0 ]]; then
    echo "[FAIL] Incomplete temporary integrated FASTQ count: ${tmp_count}" >&2
    find "$integrated_dir" -maxdepth 1 -type f -name '*.fastq.gz.tmp.*' -printf "[FAIL] %p\n" >&2
    return 1
  fi

  local audit_status
  audit_status=0
  # This is the same parallel expected-file audit used before tar packaging.
  audit_integrated_outputs "check" || audit_status=1

  echo "[INFO] Tar count: $(find "$upload_dir" -maxdepth 1 -type f -name '*.tar.gz' | wc -l)"
  if [[ -d "${upload_dir}/plain_pe" ]]; then
    echo "[INFO] Plain PE FASTQ count: $(find "${upload_dir}/plain_pe" -maxdepth 1 -type f -name '*.fastq.gz' | wc -l)"
  fi

  local tar_file fastq_count md5_count status
  status=0
  for tar_file in "$upload_dir"/*.tar.gz; do
    [[ -e "$tar_file" ]] || continue
    fastq_count="$(tar -tzf "$tar_file" | grep -c '\.fastq\.gz$' || true)"
    md5_count="$(tar -tzf "$tar_file" | grep -c '_md5\.txt$' || true)"
    if [[ "$md5_count" -ne 1 ]]; then
      echo "[FAIL] $(basename "$tar_file"): md5 file count=${md5_count}" >&2
      status=1
    elif [[ "$fastq_count" -gt 4 ]]; then
      echo "[FAIL] $(basename "$tar_file"): more than 4 FASTQs after lane integration (${fastq_count})" >&2
      status=1
    else
      echo "[OK] $(basename "$tar_file"): FASTQs=${fastq_count}, md5=${md5_count}"
    fi
  done
  if [[ "$audit_status" -ne 0 ]]; then
    status=1
  fi
  return "$status"
}

send_email_notification() {
  # Notification failure should not hide the real status of FTP/integration.
  # Missing msmtp/mail therefore warns and returns success from this helper.
  local status="$1"
  local step="$2"
  local detail="$3"
  local subject

  [[ -n "$notify_email" ]] || return 0
  if ! command -v "$mail_cmd" >/dev/null 2>&1; then
    echo "[WARN] Mail command not found, skipping email notification: ${mail_cmd}" >&2
    return 0
  fi

  subject="[ftp.sh] ${step} ${status}: $(basename "$PWD")"
  {
    printf "From: %s\n" "$mail_from"
    printf "To: %s\n" "$notify_email"
    printf "Subject: %s\n" "$subject"
    printf "\n"
    printf "Status: %s\n" "$status"
    printf "Step: %s\n" "$step"
    printf "Time: %s\n" "$(date '+%F %T')"
    printf "Host: %s\n" "$(hostname)"
    printf "Workdir: %s\n" "$PWD"
    printf "Upload directory: %s\n" "$upload_dir"
    printf "Remote: %s/%s\n" "$ftp_host" "$remote_dir"
    printf "\n%s\n" "$detail"
  } | "$mail_cmd" "$notify_email" || echo "[WARN] Failed to send email notification to ${notify_email}" >&2
}

run_step() {
  # Wrap each top-level command so success/failure emails are sent consistently.
  local step="$1"
  shift

  if "$@"; then
    send_email_notification "SUCCESS" "$step" "${step} step finished successfully."
  else
    local exit_code="$?"
    send_email_notification "FAILED" "$step" "${step} step failed with exit code ${exit_code}. Check the ftp.sh terminal output and rerun after fixing the issue."
    return "$exit_code"
  fi
}

test_email() {
  # Lightweight configuration check for msmtp/mail without running any FASTQ
  # processing or FTP transfer.
  [[ -n "$notify_email" ]] || die "Notification email is empty. Set NOTIFY_EMAIL or --email."
  send_email_notification "TEST" "test-email" "This is a test email from ftp.sh. Email notification is working."
  log "Sent test email to ${notify_email} from ${mail_from} using ${mail_cmd}."
}

upload_all() {
  # Upload only the prepared upload directory. Integration and tar are separate
  # commands, so failed uploads can be retried without touching FASTQs.
  require_cmd lftp
  [[ -n "$ftp_user" ]] || die "FTP user is empty. Set FTP_USER or --ftp-user."
  [[ -n "$ftp_pass" ]] || die "FTP password is empty. Set FTP_PASS or --ftp-pass."
  [[ -n "$remote_dir" ]] || die "Remote directory is empty. Set REMOTE_DIR or --remote-dir."
  [[ -d "$upload_dir" ]] || die "Upload directory not found: $upload_dir"

  log "Uploading ${upload_dir} to ${ftp_host}/${remote_dir}."
  lftp -u "$ftp_user","$ftp_pass" "$ftp_host" <<EOF
set ftp:ssl-allow no
set net:max-retries 3
set net:timeout 30
mkdir -p "$remote_dir"
cd "$remote_dir"
mirror -R --verbose --only-newer "$upload_dir" .
bye
EOF
}

command_name="${1:-all}"
if [[ "$command_name" == "-h" || "$command_name" == "--help" ]]; then
  usage
  exit 0
fi
shift || true
parse_args "$@"

case "$command_name" in
  integrate)
    register_running_pid
    run_step "integrate" integrate_all
    ;;
  tar)
    register_running_pid
    run_step "tar" tar_all
    ;;
  check)
    run_step "check" check_all
    ;;
  test-email)
    test_email
    ;;
  ftp|upload)
    register_running_pid
    run_step "ftp" upload_all
    ;;
  all)
    register_running_pid
    run_step "integrate" integrate_all
    run_step "tar" tar_all
    run_step "check" check_all
    run_step "ftp" upload_all
    ;;
  stop)
    stop_running_job
    ;;
  *)
    usage >&2
    die "Unknown command: $command_name"
    ;;
esac
