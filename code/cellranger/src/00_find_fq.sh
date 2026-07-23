#!/bin/bash
# ============================================================================
# 00_find_fq.sh
# Find all *.fastq.gz or *.fq.gz files under an input directory and create
# symlinks in an output directory.
#
# Usage:
#   bash cellranger/src/00_find_fq.sh --dir INPUT_DIR --outdir OUTPUT_DIR [--suffix EXT]
#
# Options:
#   --dir DIR      Root directory to search recursively for FASTQ files
#   --outdir DIR   Directory where symlinks will be created
#   --suffix EXT   Optional suffix filter: fq.gz or fastq.gz
#                  Can be passed multiple times. Default: both.
#   -h, --help     Show this help
# ============================================================================

set -euo pipefail

SEARCH_DIR=""
OUTDIR=""
SUFFIXES=("fastq.gz" "fq.gz")
CUSTOM_SUFFIX=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dir)
            SEARCH_DIR="${2:-}"
            shift 2
            ;;
        --outdir)
            OUTDIR="${2:-}"
            shift 2
            ;;
        --suffix)
            suffix="${2:-}"
            if [[ -z "${suffix}" ]]; then
                echo "[ERROR] --suffix requires a value: fq.gz or fastq.gz"
                exit 1
            fi
            if [[ "${CUSTOM_SUFFIX}" == false ]]; then
                SUFFIXES=()
                CUSTOM_SUFFIX=true
            fi
            case "${suffix}" in
                fq.gz|fastq.gz)
                    if [[ ! " ${SUFFIXES[*]} " =~ [[:space:]]${suffix}[[:space:]] ]]; then
                        SUFFIXES+=("${suffix}")
                    fi
                    ;;
                *)
                    echo "[ERROR] Invalid --suffix '${suffix}'. Allowed: fq.gz, fastq.gz"
                    exit 1
                    ;;
            esac
            shift 2
            ;;
        -h|--help)
            sed -n '/^# Usage/,/^# ====/p' "$0" | sed 's/^# \?//'
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            exit 1
            ;;
    esac
done

if [[ -z "${SEARCH_DIR}" || -z "${OUTDIR}" ]]; then
    echo "[ERROR] Both --dir and --outdir are required."
    echo "        Example: bash cellranger/src/00_find_fq.sh --dir raw_data --outdir fastq_links"
    exit 1
fi

if [[ ! -d "${SEARCH_DIR}" ]]; then
    echo "[ERROR] Input directory not found: ${SEARCH_DIR}"
    exit 1
fi

mkdir -p "${OUTDIR}"

SEARCH_DIR_ABS="$(cd "${SEARCH_DIR}" && pwd)"
OUTDIR_ABS="$(cd "${OUTDIR}" && pwd)"

found=0
linked=0
skipped=0
errors=0
fq_seen=0
source_renamed=0
rename_decided=false
rename_fq=false
clean_match_count=0
clean_removed=0

find_args=("${SEARCH_DIR_ABS}" -type f "(")
for i in "${!SUFFIXES[@]}"; do
    if [[ "${i}" -gt 0 ]]; then
        find_args+=(-o)
    fi
    find_args+=(-name "*.${SUFFIXES[$i]}")
done
find_args+=(")" -print0)

while IFS= read -r -d '' fq; do
    found=$((found + 1))

    base_name="$(basename "${fq}")"

    if [[ "${base_name}" == *.fq.gz ]]; then
        fq_seen=$((fq_seen + 1))

        if [[ "${rename_decided}" == false ]]; then
            rename_decided=true
            if [[ -t 0 ]]; then
                echo "[INFO]  Detected .fq.gz files."
                echo "[INFO]  Cell Ranger only recognizes .fastq.gz reliably."
                read -r -p "Rename source .fq.gz files to .fastq.gz? [y/N]: " reply
                case "${reply}" in
                    y|Y|yes|YES) rename_fq=true ;;
                    *) rename_fq=false ;;
                esac
            else
                echo "[WARN]  .fq.gz files detected but no interactive terminal; keeping original names."
                rename_fq=false
            fi
        fi

        if [[ "${rename_fq}" == true ]]; then
            renamed_src="${fq%.fq.gz}.fastq.gz"
            if [[ -e "${renamed_src}" || -L "${renamed_src}" ]]; then
                echo "[WARN]  Rename target exists, skipping source rename and link: ${renamed_src} (source: ${fq})"
                skipped=$((skipped + 1))
                continue
            fi
            if mv -- "${fq}" "${renamed_src}"; then
                fq="${renamed_src}"
                base_name="$(basename "${fq}")"
                source_renamed=$((source_renamed + 1))
            else
                echo "[ERROR] Failed to rename source file: ${fq} -> ${renamed_src}"
                errors=$((errors + 1))
                continue
            fi
        fi
    fi

    link_path="${OUTDIR_ABS}/${base_name}"

    if [[ -e "${link_path}" || -L "${link_path}" ]]; then
        # Keep existing correct links; avoid overwriting any existing file.
        if [[ -L "${link_path}" ]]; then
            current_target="$(readlink "${link_path}")"
            target_abs="$(realpath "${fq}")"
            current_abs="$(realpath -m "${OUTDIR_ABS}/${current_target}")"
            if [[ "${current_abs}" == "${target_abs}" ]]; then
                skipped=$((skipped + 1))
                continue
            fi
        fi
        echo "[WARN]  Exists, skipping: ${link_path}"
        skipped=$((skipped + 1))
        continue
    fi

    if ln -s "${fq}" "${link_path}"; then
        linked=$((linked + 1))
    else
        echo "[ERROR] Failed to create symlink: ${link_path} -> ${fq}"
        errors=$((errors + 1))
    fi
done < <(find "${find_args[@]}")

# Post-link hygiene check: optionally remove symlinks whose link path or target
# contains Clean/clean because these are often undesired duplicate FASTQ sets.
clean_links=()
keep_links=()
while IFS= read -r -d '' link_path; do
    link_target="$(readlink "${link_path}")"
    target_abs="$(realpath -m "${OUTDIR_ABS}/${link_target}")"
    if [[ "${link_path}" == *[Cc]lean* || "${link_target}" == *[Cc]lean* || "${target_abs}" == *[Cc]lean* ]]; then
        clean_links+=("${link_path}")
    else
        keep_links+=("${link_path}")
    fi
done < <(find "${OUTDIR_ABS}" -maxdepth 1 -type l -print0)

clean_match_count=${#clean_links[@]}
if [[ "${clean_match_count}" -gt 0 ]]; then
    echo "[WARN]  Found ${clean_match_count} symlink(s) with 'Clean/clean' in link path or target."
    echo "[INFO]  Links suggested for removal:"
    for lnk in "${clean_links[@]}"; do
        echo "        ${lnk}"
    done

    echo "[INFO]  Links to keep:"
    if [[ "${#keep_links[@]}" -eq 0 ]]; then
        echo "        (none)"
    else
        for lnk in "${keep_links[@]}"; do
            echo "        ${lnk}"
        done
    fi

    if [[ -t 0 ]]; then
        read -r -p "Remove the listed Clean/clean symlinks from outdir? [y/N]: " clean_reply
        case "${clean_reply}" in
            y|Y|yes|YES)
                for lnk in "${clean_links[@]}"; do
                    if rm -- "${lnk}"; then
                        clean_removed=$((clean_removed + 1))
                    else
                        echo "[ERROR] Failed to remove symlink: ${lnk}"
                        errors=$((errors + 1))
                    fi
                done
                ;;
            *)
                echo "[INFO]  Keeping listed Clean/clean symlinks."
                ;;
        esac
    else
        echo "[WARN]  Non-interactive mode; keeping listed Clean/clean symlinks."
    fi
fi

echo "[INFO]  Search dir:      ${SEARCH_DIR_ABS}"
echo "[INFO]  Out dir:         ${OUTDIR_ABS}"
echo "[INFO]  Suffixes:        ${SUFFIXES[*]}"
echo "[INFO]  Found:           ${found}"
echo "[INFO]  Linked:          ${linked}"
echo "[INFO]  Skipped:         ${skipped}"
echo "[INFO]  Errors:          ${errors}"
echo "[INFO]  fq.gz seen:      ${fq_seen}"
echo "[INFO]  Source renamed:  ${source_renamed}"
echo "[INFO]  Rename enabled:  ${rename_fq}"
echo "[INFO]  Clean matches:   ${clean_match_count}"
echo "[INFO]  Clean removed:   ${clean_removed}"
echo "[INFO]  Done."

if [[ "${found}" -eq 0 ]]; then
    echo "[WARN]  No matching FASTQ files found under ${SEARCH_DIR_ABS}"
fi

if [[ "${errors}" -gt 0 ]]; then
    exit 1
fi
