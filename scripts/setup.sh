#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PARAMS="${ROOT_DIR}/code/config/params.yaml"
DRY_RUN="false"
SKIP_EXISTING="true"
PARTS=()

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/setup.sh [all|env|software|ref|r|check] [options]

Options:
  --params PATH       params.yaml path
  --part NAME         Setup one part; may be repeated. Names: env, software, ref, r
  --dry-run           Print commands without running them
  --no-skip-existing  Re-extract local archives even when target paths exist
  -h, --help          Show help

Examples:
  bash scripts/setup.sh check
  bash scripts/setup.sh env --params code/config/params.yaml
  bash scripts/setup.sh --part software --part ref
  bash scripts/setup.sh all --dry-run
USAGE
}

log() {
  echo "[$(date '+%F %T')] $*"
}

die() {
  echo "[ERROR] $*" >&2
  exit 1
}

run_cmd() {
  if [[ "$DRY_RUN" == "true" ]]; then
    printf '[DRY-RUN]'
    printf ' %q' "$@"
    printf '\n'
  else
    "$@"
  fi
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

yaml_get() {
  local key="$1"
  local default="${2:-}"
  python3 - "$PARAMS" "$key" "$default" <<'PY'
import sys

params, key, default = sys.argv[1:4]
try:
    import yaml
except Exception as exc:
    print(f"[ERROR] Python package PyYAML is required to read {params}: {exc}", file=sys.stderr)
    sys.exit(2)

with open(params, encoding="utf-8") as handle:
    data = yaml.safe_load(handle) or {}

value = data
for part in key.split("."):
    if not isinstance(value, dict) or part not in value:
        print(default)
        sys.exit(0)
    value = value[part]

if value is None:
    print(default)
elif isinstance(value, bool):
    print("true" if value else "false")
elif isinstance(value, (list, tuple)):
    print("\n".join(str(x) for x in value))
else:
    print(str(value))
PY
}

yaml_list() {
  yaml_get "$1" "" | awk 'NF'
}

is_abs_path() {
  [[ "${1:-}" == /* ]]
}

abs_or_root_path() {
  local path="$1"
  if [[ -z "$path" ]]; then
    echo ""
  elif is_abs_path "$path"; then
    echo "$path"
  else
    echo "${ROOT_DIR}/${path}"
  fi
}

ensure_dir() {
  local path="$1"
  [[ -n "$path" ]] || return 0
  if [[ -d "$path" ]]; then
    log "Directory exists: $path"
  else
    run_cmd mkdir -p "$path"
  fi
}

extract_archive() {
  local archive="$1"
  local dest="$2"
  [[ -n "$archive" ]] || return 1
  [[ -f "$archive" ]] || die "Archive not found: $archive"
  ensure_dir "$dest"
  case "$archive" in
    *.tar.gz|*.tgz) run_cmd tar -xzf "$archive" -C "$dest" ;;
    *.tar.bz2|*.tbz2) run_cmd tar -xjf "$archive" -C "$dest" ;;
    *.tar.xz|*.txz) run_cmd tar -xJf "$archive" -C "$dest" ;;
    *.tar) run_cmd tar -xf "$archive" -C "$dest" ;;
    *.zip)
      require_cmd unzip
      run_cmd unzip -q "$archive" -d "$dest"
      ;;
    *) die "Unsupported archive type: $archive" ;;
  esac
}

link_executable() {
  local source="$1"
  local link="$2"
  [[ -x "$source" ]] || die "Executable not found or not executable: $source"
  ensure_dir "$(dirname "$link")"
  run_cmd ln -sfn "$source" "$link"
}

find_first_executable() {
  local root="$1"
  local name="$2"
  find "$root" -type f -name "$name" -perm /111 2>/dev/null | sort | head -n 1
}

setup_env() {
  local outdir project_dir bin_dir r_lib
  outdir="$(yaml_get outdir output)"
  project_dir="$(yaml_get setup.project_dir "$ROOT_DIR")"
  bin_dir="$(abs_or_root_path "$(yaml_get setup.bin_dir tools/bin)")"
  r_lib="$(abs_or_root_path "$(yaml_get setup.r_lib renv/library)")"

  log "Setting up project directories."
  ensure_dir "$project_dir"
  ensure_dir "$outdir"
  ensure_dir "$bin_dir"
  ensure_dir "$r_lib"
  ensure_dir "${ROOT_DIR}/rds"
  ensure_dir "${ROOT_DIR}/figures"
  ensure_dir "${ROOT_DIR}/table"
  ensure_dir "${ROOT_DIR}/log"
  ensure_dir "${ROOT_DIR}/tmp"

  local conda_prefix conda_yaml
  conda_prefix="$(yaml_get setup.conda_env_prefix "")"
  conda_yaml="$(yaml_get setup.conda_env_yaml "")"
  if [[ -n "$conda_prefix" && -n "$conda_yaml" ]]; then
    require_cmd conda
    if [[ -d "$conda_prefix" && "$SKIP_EXISTING" == "true" ]]; then
      log "Conda env exists: $conda_prefix"
    else
      run_cmd conda env create -p "$conda_prefix" -f "$conda_yaml"
    fi
  fi
}

setup_software() {
  local software_dir bin_dir cellranger_path fastqc_path cellranger_archive fastqc_archive
  software_dir="$(yaml_get setup.software_dir /path/to/software)"
  bin_dir="$(abs_or_root_path "$(yaml_get setup.bin_dir tools/bin)")"
  cellranger_path="$(yaml_get cellranger.cellranger_path "")"
  fastqc_path="$(yaml_get cellranger.fastqc_path "")"
  cellranger_archive="$(yaml_get setup.cellranger_archive "")"
  fastqc_archive="$(yaml_get setup.fastqc_archive "")"

  ensure_dir "$software_dir"
  ensure_dir "$bin_dir"

  if [[ -x "$cellranger_path" && "$SKIP_EXISTING" == "true" ]]; then
    log "Cell Ranger exists: $cellranger_path"
  elif [[ -n "$cellranger_archive" ]]; then
    extract_archive "$cellranger_archive" "$software_dir"
    cellranger_path="$(find_first_executable "$software_dir" cellranger)"
  fi
  [[ -x "$cellranger_path" ]] && link_executable "$cellranger_path" "${bin_dir}/cellranger" || log "Cell Ranger not installed; set cellranger.cellranger_path or setup.cellranger_archive."

  if [[ -x "$fastqc_path" && "$SKIP_EXISTING" == "true" ]]; then
    log "FastQC exists: $fastqc_path"
  elif [[ -n "$fastqc_archive" ]]; then
    extract_archive "$fastqc_archive" "$software_dir"
    fastqc_path="$(find_first_executable "$software_dir" fastqc)"
  fi
  [[ -x "$fastqc_path" ]] && link_executable "$fastqc_path" "${bin_dir}/fastqc" || log "FastQC not installed; set cellranger.fastqc_path or setup.fastqc_archive."
}

setup_ref() {
  local reference_dir gex_ref vdj_ref gex_archive vdj_archive
  reference_dir="$(yaml_get setup.reference_dir /path/to/reference/cellranger)"
  gex_ref="$(yaml_get cellranger.reference_gex "")"
  vdj_ref="$(yaml_get cellranger.reference_vdj "")"
  gex_archive="$(yaml_get setup.gex_reference_archive "")"
  vdj_archive="$(yaml_get setup.vdj_reference_archive "")"

  ensure_dir "$reference_dir"
  if [[ -d "$gex_ref" && "$SKIP_EXISTING" == "true" ]]; then
    log "GEX reference exists: $gex_ref"
  elif [[ -n "$gex_archive" ]]; then
    extract_archive "$gex_archive" "$reference_dir"
  else
    log "GEX reference not installed; set cellranger.reference_gex or setup.gex_reference_archive."
  fi

  if [[ -d "$vdj_ref" && "$SKIP_EXISTING" == "true" ]]; then
    log "VDJ reference exists: $vdj_ref"
  elif [[ -n "$vdj_archive" ]]; then
    extract_archive "$vdj_archive" "$reference_dir"
  else
    log "VDJ reference not installed; set cellranger.reference_vdj or setup.vdj_reference_archive."
  fi
}

setup_r() {
  require_cmd Rscript
  local r_lib package_source lockfile
  r_lib="$(abs_or_root_path "$(yaml_get setup.r_lib renv/library)")"
  package_source="$(yaml_get setup.r_package_source_dir "")"
  lockfile="${ROOT_DIR}/renv.lock"
  ensure_dir "$r_lib"

  if [[ "$DRY_RUN" == "true" ]]; then
    log "Would check/install R packages from ${lockfile} into ${r_lib}; local package source: ${package_source:-<unset>}."
    return 0
  fi

  Rscript --vanilla - "$lockfile" "$r_lib" "$package_source" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
lockfile <- args[[1]]
lib <- args[[2]]
source_dir <- args[[3]]

dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

lock <- jsonlite::fromJSON(lockfile)
packages <- names(lock$Packages)
packages <- setdiff(packages, c("base", "stats", "grid", "parallel"))
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

if (!length(missing)) {
  message("[OK] Required R packages already available in .libPaths().")
  quit(status = 0)
}

if (!nzchar(source_dir) || !dir.exists(source_dir)) {
  message("[WARN] Missing R packages: ", paste(missing, collapse = ", "))
  message("[WARN] setup.r_package_source_dir is empty or missing; local install skipped.")
  quit(status = 0)
}

archives <- list.files(source_dir, pattern = "\\.(tar\\.gz|tgz|zip)$", full.names = TRUE)
for (pkg in missing) {
  hit <- archives[grepl(paste0("^", pkg, "_"), basename(archives))]
  if (!length(hit)) {
    message("[WARN] No local archive found for R package: ", pkg)
    next
  }
  message("[RUN] Installing R package from local archive: ", basename(hit[[1]]))
  install.packages(hit[[1]], repos = NULL, type = "source", lib = lib)
}
RS
}

check_setup() {
  log "Checking configured paths."
  local path
  for path in \
    "$(yaml_get cellranger.cellranger_path "")" \
    "$(yaml_get cellranger.fastqc_path "")" \
    "$(yaml_get cellranger.reference_gex "")" \
    "$(yaml_get cellranger.reference_vdj "")" \
    "$(yaml_get cellranger.inner_enrichment_primers "")"; do
    [[ -n "$path" ]] || continue
    if [[ -e "$path" ]]; then
      echo "[OK] $path"
    else
      echo "[MISS] $path"
    fi
  done
}

run_part() {
  case "$1" in
    env) setup_env ;;
    software) setup_software ;;
    ref) setup_ref ;;
    r) setup_r ;;
    check) check_setup ;;
    *) die "Unknown setup part: $1" ;;
  esac
}

command_name="${1:-all}"
if [[ "$command_name" == "-h" || "$command_name" == "--help" ]]; then
  usage
  exit 0
fi
shift || true

while [[ $# -gt 0 ]]; do
  case "$1" in
    --params) PARAMS="$2"; shift 2 ;;
    --part) PARTS+=("$2"); shift 2 ;;
    --dry-run) DRY_RUN="true"; shift ;;
    --no-skip-existing) SKIP_EXISTING="false"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -f "$PARAMS" ]] || die "params.yaml not found: $PARAMS"
require_cmd python3

case "$command_name" in
  all)
    if [[ "${#PARTS[@]}" -eq 0 ]]; then
      mapfile -t PARTS < <(yaml_list setup.parts)
      [[ "${#PARTS[@]}" -gt 0 ]] || PARTS=(env software ref r)
    fi
    ;;
  env|software|ref|r|check)
    PARTS=("$command_name")
    ;;
  *)
    die "Unknown command: $command_name"
    ;;
esac

for part in "${PARTS[@]}"; do
  run_part "$part"
done
