#!/usr/bin/env python3
"""Generate per-sample Cell Ranger multi config CSVs."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys
from pathlib import Path
from typing import Any

try:
    import yaml
except Exception as exc:  # pragma: no cover - dependency validation path
    print(f"[ERROR] Python package 'yaml' (PyYAML) is required: {exc}", file=sys.stderr)
    sys.exit(1)


DEFAULT_PARAMS = "config/params.yaml"
DEFAULT_SAMPLES = "config/samples.csv"
DEFAULT_TEMPLATE = "cellranger/src/multi_template.csv"
DEFAULT_OUTDIR = "output"


def find_project_root(start: Path) -> Path:
    current = start.resolve()
    while True:
        if (current / "config").is_dir() and (current / "cellranger" / "src").is_dir():
            return current
        if current.parent == current:
            raise SystemExit(f"[ERROR] Could not locate project root from: {start}")
        current = current.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate per-sample Cell Ranger multi config CSVs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python3 cellranger/src/02_generate_cellranger_configs.py \\\n"
            "    --params config/params.yaml \\\n"
            "    --samples config/samples.csv \\\n"
            "    --template cellranger/src/multi_template.csv \\\n"
            "    --outdir output"
        ),
    )
    parser.add_argument("--params", default=DEFAULT_PARAMS, help=f"Path to params.yaml [default: {DEFAULT_PARAMS}]")
    parser.add_argument("--samples", default=DEFAULT_SAMPLES, help=f"Path to samples.csv [default: {DEFAULT_SAMPLES}]")
    parser.add_argument("--template", default=DEFAULT_TEMPLATE, help=f"Path to multi_template.csv [default: {DEFAULT_TEMPLATE}]")
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR, help=f"Pipeline output root [default: {DEFAULT_OUTDIR}]")
    parser.add_argument("--dry-run", action="store_true", help="Use defaults and skip path validation")
    return parser.parse_args()


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise SystemExit(f"[ERROR] Top-level YAML must be a mapping: {path}")
    return data


def load_samples(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        missing = [col for col in ("sample_id", "patient_id", "group") if col not in fields]
        if missing:
            raise SystemExit("[ERROR] samples.csv missing required columns: " + ", ".join(missing))

        rows: list[dict[str, str]] = []
        seen: set[str] = set()
        for line_no, raw_row in enumerate(reader, start=2):
            row = {key: (value or "").strip() for key, value in raw_row.items()}
            sample_id = row.get("sample_id", "")
            if not sample_id:
                raise SystemExit(f"[ERROR] samples.csv line {line_no}: sample_id is empty")
            if sample_id in seen:
                raise SystemExit(f"[ERROR] samples.csv line {line_no}: duplicate sample_id {sample_id}")
            seen.add(sample_id)
            rows.append(row)

    if not rows:
        raise SystemExit(f"[ERROR] samples.csv contains no samples: {path}")
    return rows


def is_placeholder(value: str | None) -> bool:
    return bool(value) and value.strip().startswith("/path/to/")


def is_present(value: str | None) -> bool:
    return bool(value and value.strip())


def validate_inputs(config: dict[str, Any], samples: list[dict[str, str]]) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    cellranger = config.get("cellranger") or {}
    if not isinstance(cellranger, dict):
        errors.append("Missing or invalid params key: cellranger")
        cellranger = {}

    def check_dir(path: str | None, label: str) -> None:
        if not is_present(path) or is_placeholder(path):
            errors.append(f"{label} is not set")
        elif not Path(path).is_dir():
            errors.append(f"{label} not found: {path}")

    def check_file(path: str | None, label: str, required: bool = True) -> None:
        if not is_present(path):
            if required:
                errors.append(f"{label} is not set")
            return
        if is_placeholder(path):
            errors.append(f"{label} is not set")
        elif not Path(path).is_file():
            errors.append(f"{label} not found: {path}")

    def check_fastq_dir(path: str | None, sid: str, label: str) -> None:
        if not is_present(path):
            return
        if is_placeholder(path):
            errors.append(f"samples.csv [{sid}]: {label} is not set")
        elif not Path(path).is_dir():
            errors.append(f"samples.csv [{sid}]: {label} not found: {path}")

    check_dir(str(cellranger.get("reference_gex", "") or ""), "cellranger.reference_gex")

    vdj_path_cols = ("vdj_t_fastq_path", "vdj_b_fastq_path", "vdj_t_gd_fastq_path")
    any_vdj = any(any(is_present(row.get(col)) for col in vdj_path_cols) for row in samples)
    if any_vdj:
        check_dir(str(cellranger.get("reference_vdj", "") or ""), "cellranger.reference_vdj")

    has_gd = any(is_present(row.get("vdj_t_gd_fastq_path")) for row in samples)
    if has_gd:
        check_file(
            str(cellranger.get("inner_enrichment_primers", "") or ""),
            "cellranger.inner_enrichment_primers",
            required=True,
        )

    cr_path = str(cellranger.get("cellranger_path", "") or "").strip()
    if cr_path and not is_placeholder(cr_path) and cr_path != "cellranger":
        if os.path.sep in cr_path:
            if not Path(cr_path).exists():
                warnings.append(f"cellranger.cellranger_path not found: {cr_path}")
        elif shutil.which(cr_path) is None:
            warnings.append(f"cellranger command not found on PATH: {cr_path}")

    for row in samples:
        sid = row["sample_id"]
        check_fastq_dir(row.get("gex_fastq_path"), sid, "gex_fastq_path")
        check_fastq_dir(row.get("vdj_t_fastq_path"), sid, "vdj_t_fastq_path")
        check_fastq_dir(row.get("vdj_b_fastq_path"), sid, "vdj_b_fastq_path")
        check_fastq_dir(row.get("vdj_t_gd_fastq_path"), sid, "vdj_t_gd_fastq_path")

    return errors, warnings


def resolve_pipeline_outdirs(config: dict[str, Any], outdir: str) -> dict[str, Path]:
    root = Path(outdir or config.get("outdir") or "output")
    return {
        "root": root,
        "cellranger_configs": root / "cellranger_configs",
        "cellranger_output": root / "cellranger_output",
    }


def clean(value: Any) -> str:
    return "" if value is None else str(value)


def without_vdj_block(lines: list[str]) -> list[str]:
    out: list[str] = []
    in_vdj = False
    for line in lines:
        if line == "[vdj]":
            in_vdj = True
            continue
        if in_vdj and line.startswith("["):
            in_vdj = False
        if not in_vdj:
            out.append(line)
    return out


def generate_configs(
    config: dict[str, Any],
    samples: list[dict[str, str]],
    template_path: Path,
    output_dir: Path,
    cellranger_output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    cellranger_output_dir.mkdir(parents=True, exist_ok=True)

    template = template_path.read_text(encoding="utf-8").splitlines()
    cellranger = config.get("cellranger") or {}
    replacements_common = {
        "__GEX_REF__": clean(cellranger.get("reference_gex")),
        "__VDJ_REF__": clean(cellranger.get("reference_vdj")),
        "__VDJ_INNER_ENRICHMENT_PRIMERS__": clean(cellranger.get("inner_enrichment_primers")),
    }

    for row in samples:
        replacements = {
            **replacements_common,
            "__GEX_FASTQ_PREFIX__": clean(row.get("gex_fastq_prefix")),
            "__GEX_FASTQ_PATH__": clean(row.get("gex_fastq_path")),
            "__VDJ_T_FASTQ_PREFIX__": clean(row.get("vdj_t_fastq_prefix")),
            "__VDJ_T_FASTQ_PATH__": clean(row.get("vdj_t_fastq_path")),
            "__VDJ_B_FASTQ_PREFIX__": clean(row.get("vdj_b_fastq_prefix")),
            "__VDJ_B_FASTQ_PATH__": clean(row.get("vdj_b_fastq_path")),
            "__VDJ_T_GD_FASTQ_PREFIX__": clean(row.get("vdj_t_gd_fastq_prefix")),
            "__VDJ_T_GD_FASTQ_PATH__": clean(row.get("vdj_t_gd_fastq_path")),
        }

        lines = template[:]
        for token, value in replacements.items():
            lines = [line.replace(token, value) for line in lines]

        has_vdj_t = is_present(row.get("vdj_t_fastq_path"))
        has_vdj_b = is_present(row.get("vdj_b_fastq_path"))
        has_vdj_gd = is_present(row.get("vdj_t_gd_fastq_path"))

        if not has_vdj_t:
            lines = [line for line in lines if not line.endswith(",VDJ-T,")]
        if not has_vdj_b:
            lines = [line for line in lines if not line.endswith(",VDJ-B,")]
        if not has_vdj_gd:
            lines = [line for line in lines if not line.endswith(",VDJ-T-GD,")]

        if not has_vdj_t and not has_vdj_b and not has_vdj_gd:
            lines = without_vdj_block(lines)
            print(f"  [INFO] {row['sample_id']}: no VDJ paths, GEX-only config")
        elif not has_vdj_gd:
            lines = [line for line in lines if not line.startswith("inner-enrichment-primers,")]

        out_path = output_dir / f"{row['sample_id']}_multi_config.csv"
        out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"Generated: {out_path}")


def main() -> int:
    args = parse_args()
    project_root = find_project_root(Path(__file__).resolve().parent)

    if args.dry_run:
        print("[TEST]  --dry-run: using defaults, skipping validation.")
        args.params = DEFAULT_PARAMS
        args.samples = DEFAULT_SAMPLES
        args.template = DEFAULT_TEMPLATE
        args.outdir = DEFAULT_OUTDIR
        print(f"[TEST]  params   = {args.params}")
        print(f"[TEST]  samples  = {args.samples}")
        print(f"[TEST]  template = {args.template}")
        print(f"[TEST]  outdir   = {args.outdir}")

    params_path = Path(args.params).resolve()
    samples_path = Path(args.samples).resolve()
    template_path = Path(args.template).resolve()
    os.chdir(project_root)
    for path, label in (
        (params_path, "params YAML"),
        (samples_path, "samples CSV"),
        (template_path, "multi template CSV"),
    ):
        if not path.is_file():
            print(f"[ERROR] {label} not found: {path}", file=sys.stderr)
            return 1

    config = load_yaml(params_path)
    samples = load_samples(samples_path)

    if args.dry_run:
        print("[INFO]  --dry-run: skipping path validation.")
    else:
        errors, warnings = validate_inputs(config, samples)
        if warnings:
            print("WARNINGS:")
            for warning in warnings:
                print(f"  [WARN] {warning}")
        if errors:
            print("ERRORS:")
            for error in errors:
                print(f"  [ERROR] {error}")
            print(f"[ERROR] {len(errors)} validation error(s); fix the above before proceeding.", file=sys.stderr)
            return 1
        print("Validation passed.")

    dirs = resolve_pipeline_outdirs(config, args.outdir)
    generate_configs(
        config=config,
        samples=samples,
        template_path=template_path,
        output_dir=dirs["cellranger_configs"],
        cellranger_output_dir=dirs["cellranger_output"],
    )
    print(f"Done. Generated {len(samples)} config(s) in {dirs['cellranger_configs']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
