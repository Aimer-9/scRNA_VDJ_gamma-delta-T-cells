#!/usr/bin/env python3
"""Run and validate the Python portion of the project pySCENIC workflow.

The shell entry point, ``11_pySCENIC.sh``, performs balanced cell selection
and exports a cells-by-genes CSV from Seurat. This module then:

1. validates the pinned Python environment and TF/gene overlap;
2. runs GRN inference on the CSV expression matrix;
3. runs cisTarget motif pruning on the GRN adjacencies;
4. converts the expression CSV to a chunked loom because this pySCENIC
   version requires loom input for AUCell;
5. runs AUCell and verifies that its output is a valid HDF5 loom containing
   CellID and RegulonsAUC.

Each pySCENIC CLI stage writes to a separate log when --log-dir is supplied.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run pySCENIC GRN, ctx, and AUCell steps.")
    parser.add_argument("--expr-csv", required=True)
    parser.add_argument(
        "--expr-loom",
        default="",
        help=(
            "Expression loom used by AUCell. If omitted, expression_matrix.loom "
            "is created beside --auc-loom."
        ),
    )
    parser.add_argument(
        "--loom-chunk-cells",
        type=int,
        default=250,
        help="Cells per chunk when converting expression CSV to loom.",
    )
    parser.add_argument("--tf-list", required=True)
    parser.add_argument("--ranking-db", action="append", default=[], required=True)
    parser.add_argument("--motif-annotations", required=True)
    parser.add_argument("--adj-tsv", required=True)
    parser.add_argument("--regulons-csv", required=True)
    parser.add_argument("--auc-loom", required=True)
    parser.add_argument("--num-workers", default="16")
    parser.add_argument(
        "--grn-num-workers",
        default="",
        help=(
            "Worker count for pySCENIC grn. If omitted, defaults to min(--num-workers, 4); "
            "the shell wrapper defaults to 1 for robust 10k-cell runs."
        ),
    )
    parser.add_argument(
        "--grn-method",
        default="grnboost2",
        choices=["grnboost2", "genie3"],
        help="pySCENIC grn inference method.",
    )
    parser.add_argument("--seed", default="", help="Optional pySCENIC grn random seed.")
    parser.add_argument(
        "--ctx-mode",
        default="custom_multiprocessing",
        choices=["custom_multiprocessing", "dask_multiprocessing"],
        help="pySCENIC ctx execution mode. custom_multiprocessing avoids Dask scheduler errors.",
    )
    parser.add_argument("--pyscenic-command", default="pyscenic")
    parser.add_argument("--log-dir", default="", help="Directory for per-step pySCENIC logs.")
    return parser.parse_args()


def check_file(path: str, label: str) -> None:
    file_path = Path(path)
    if not file_path.is_file() or file_path.stat().st_size == 0:
        raise SystemExit(f"Missing {label}: {path}")


def read_expression_genes(expr_csv: str) -> set[str]:
    with open(expr_csv, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise SystemExit(f"Expression matrix is empty: {expr_csv}") from exc
    genes = [value.strip() for value in header[1:] if value.strip()]
    if not genes:
        raise SystemExit(
            "Expression matrix header has no gene columns. "
            "pySCENIC expects cells as rows and genes as columns."
        )
    return set(genes)


def read_expression_gene_order(expr_csv: str) -> list[str]:
    with open(expr_csv, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise SystemExit(f"Expression matrix is empty: {expr_csv}") from exc
    genes = [value.strip() for value in header[1:]]
    if not genes or any(not gene for gene in genes):
        raise SystemExit(
            "Expression CSV must have a cell-ID first column followed by named gene columns."
        )
    if len(genes) != len(set(genes)):
        raise SystemExit("Expression CSV contains duplicate gene columns.")
    return genes


def validate_expression_loom(expr_loom: str) -> None:
    try:
        import h5py  # noqa: PLC0415

        with h5py.File(expr_loom, "r") as loom:
            missing = [
                path
                for path in ("matrix", "row_attrs/Gene", "col_attrs/CellID")
                if path not in loom
            ]
            shape = loom["matrix"].shape if "matrix" in loom else None
    except Exception as exc:
        raise SystemExit(f"Expression loom cannot be opened: {expr_loom}") from exc

    if missing:
        raise SystemExit(
            f"Expression loom is missing required dataset(s): {', '.join(missing)}"
        )
    if shape is None or shape[0] == 0 or shape[1] == 0:
        raise SystemExit(f"Expression loom has an empty matrix: {expr_loom}")

    print(
        f"[INFO] Expression loom matrix: {shape[0]} genes x {shape[1]} cells",
        flush=True,
    )


def convert_expression_csv_to_loom(
    expr_csv: str, expr_loom: str, chunk_cells: int
) -> None:
    """Convert cells-by-genes CSV to genes-by-cells loom without loading all cells.

    pandas reads a bounded number of cell rows per iteration. Each chunk is
    transposed to the loom convention and appended as columns. Gene order must
    remain identical across chunks because loom row attributes are written
    only with the first chunk.
    """
    if chunk_cells < 1:
        raise SystemExit("--loom-chunk-cells must be >= 1")

    try:
        import loompy  # noqa: PLC0415
        import numpy as np  # noqa: PLC0415
        import pandas as pd  # noqa: PLC0415
    except ImportError as exc:
        raise SystemExit(
            "CSV-to-loom conversion requires loompy, numpy, and pandas in the "
            "pySCENIC environment."
        ) from exc

    genes = np.asarray(read_expression_gene_order(expr_csv), dtype=object)
    loom_path = Path(expr_loom)
    loom_path.parent.mkdir(parents=True, exist_ok=True)
    if loom_path.exists():
        loom_path.unlink()

    print(
        f"[INFO] Converting expression CSV to loom in {chunk_cells}-cell chunks: "
        f"{expr_loom}",
        flush=True,
    )

    cells_written = 0
    with loompy.new(str(loom_path)) as loom:
        for chunk_index, chunk in enumerate(
            pd.read_csv(expr_csv, chunksize=chunk_cells)
        ):
            if chunk.shape[1] != len(genes) + 1:
                raise SystemExit(
                    "Expression CSV column count changed while converting to loom."
                )
            cell_ids = chunk.iloc[:, 0].astype(str).to_numpy()
            matrix = chunk.iloc[:, 1:].to_numpy(dtype=np.float32, copy=False).T
            row_attrs = {"Gene": genes} if chunk_index == 0 else None
            loom.add_columns(
                matrix,
                {"CellID": cell_ids},
                row_attrs=row_attrs,
            )
            cells_written += len(cell_ids)

    if cells_written == 0:
        raise SystemExit(f"Expression CSV has no cell rows: {expr_csv}")
    validate_expression_loom(expr_loom)


def read_tf_list(tf_list: str) -> set[str]:
    with open(tf_list, encoding="utf-8") as handle:
        tfs = {
            line.strip().split()[0]
            for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        }
    if not tfs:
        raise SystemExit(f"TF list is empty: {tf_list}")
    return tfs


def validate_tf_overlap(expr_csv: str, tf_list: str) -> None:
    genes = read_expression_genes(expr_csv)
    tfs = read_tf_list(tf_list)
    overlap = genes.intersection(tfs)
    if not overlap:
        raise SystemExit(
            "No overlap between expression-matrix gene columns and TF_LIST. "
            "This often appears later as a Dask error like "
            "'Must supply at least one delayed object'. Check gene symbols, "
            "species, and TF list file."
        )
    print(f"[INFO] TF overlap with expression genes: {len(overlap)} / {len(tfs)}", flush=True)


def count_csv_rows(path: str, delimiter: str = ",") -> int:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        try:
            next(reader)
        except StopIteration:
            return 0
        return sum(1 for _ in reader)


def validate_grn_output(adj_tsv: str, expr_csv: str, tf_list: str, log_path: Path | None) -> None:
    adj_path = Path(adj_tsv)
    if not adj_path.is_file() or adj_path.stat().st_size == 0:
        log_hint = f" Inspect GRN log: {log_path}" if log_path else ""
        raise SystemExit(
            f"pySCENIC GRN created no adjacencies: {adj_tsv}.{log_hint} "
            "Common causes are Dask worker failure, too many GRN workers, or no usable TFs "
            "after expression filtering. Retry with --grn-num-workers 1 and --seed 1."
        )

    n_edges = count_csv_rows(adj_tsv, delimiter="\t")
    if n_edges == 0:
        genes = read_expression_genes(expr_csv)
        tfs = read_tf_list(tf_list)
        overlap = genes.intersection(tfs)
        log_hint = f" Inspect GRN log: {log_path}" if log_path else ""
        raise SystemExit(
            f"pySCENIC GRN output has zero TF-target edges: {adj_tsv}.{log_hint} "
            f"Expression genes: {len(genes)}; TF list entries: {len(tfs)}; "
            f"overlapping TFs: {len(overlap)}. "
            "Check that gene symbols match the TF list species/case and that expression filtering "
            "did not remove most TF genes. If overlap is reasonable, retry with "
            "--grn-num-workers 1 --seed 1 or try --grn-method genie3."
        )

    print(f"[INFO] GRN adjacency edges: {n_edges}", flush=True)


def validate_auc_loom(auc_loom: str, log_path: Path | None) -> None:
    loom_path = Path(auc_loom)
    log_hint = f" Inspect AUCell log: {log_path}" if log_path else ""
    if not loom_path.is_file() or loom_path.stat().st_size < 8:
        raise SystemExit(
            f"pySCENIC AUCell did not create a valid loom file: {auc_loom}.{log_hint}"
        )

    with loom_path.open("rb") as handle:
        signature = handle.read(8)
    if signature != b"\x89HDF\r\n\x1a\n":
        preview = loom_path.read_bytes()[:160].decode("utf-8", errors="replace").strip()
        preview_hint = f" File begins with: {preview!r}." if preview else ""
        raise SystemExit(
            f"pySCENIC AUCell output is not an HDF5/loom file: {auc_loom}."
            f"{preview_hint}{log_hint} Remove the invalid file and rerun AUCell."
        )

    try:
        import h5py  # noqa: PLC0415

        with h5py.File(loom_path, "r") as loom:
            missing = [
                path
                for path in ("col_attrs/CellID", "col_attrs/RegulonsAUC")
                if path not in loom
            ]
    except Exception as exc:
        raise SystemExit(
            f"AUCell output has an HDF5 signature but cannot be opened: {auc_loom}."
            f"{log_hint}"
        ) from exc

    if missing:
        raise SystemExit(
            f"AUCell loom is missing required dataset(s): {', '.join(missing)}."
            f"{log_hint}"
        )

    print(
        f"[INFO] Valid AUCell loom: {auc_loom} ({loom_path.stat().st_size} bytes)",
        flush=True,
    )


def validate_pyscenic_python_environment() -> None:
    try:
        import dask  # noqa: PLC0415
        import distributed  # noqa: PLC0415
        import pandas as pd  # noqa: PLC0415
    except ImportError as exc:
        raise SystemExit(
            "Missing Python package in the pySCENIC environment. "
            "Create/activate code/downstream/envs/pyscenic_environment.yml "
            "or pass --pyscenic-python from that environment."
        ) from exc

    if not hasattr(pd.core.strings, "StringMethods"):
        raise SystemExit(
            "Incompatible pandas version for pySCENIC/Dask: "
            f"pandas {pd.__version__} from {Path(pd.__file__).parent}. "
            "This environment triggers "
            "AttributeError: module 'pandas.core.strings' has no attribute 'StringMethods'. "
            "Use the pinned pySCENIC environment with pandas==1.5.3, dask==2022.12.1, "
            "and distributed==2022.12.1."
        )

    print(
        "[INFO] Python package versions: "
        f"pandas={pd.__version__}, dask={dask.__version__}, "
        f"distributed={distributed.__version__}",
        flush=True,
    )


def parse_positive_int(value: str, label: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise SystemExit(f"{label} must be an integer: {value}") from exc
    if parsed < 1:
        raise SystemExit(f"{label} must be >= 1: {value}")
    return parsed


def default_grn_workers(num_workers: str) -> str:
    return str(min(parse_positive_int(num_workers, "--num-workers"), 4))


def configure_subprocess_threads() -> None:
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        if not os.environ.get(variable):
            os.environ[variable] = "1"

    semaphore_warning_filter = (
        "ignore:resource_tracker:UserWarning:multiprocessing.resource_tracker"
    )
    existing_warnings = os.environ.get("PYTHONWARNINGS", "")
    warning_filters = [item for item in existing_warnings.split(",") if item]
    if semaphore_warning_filter not in warning_filters:
        warning_filters.append(semaphore_warning_filter)
        os.environ["PYTHONWARNINGS"] = ",".join(warning_filters)


def resolve_pyscenic_command(command: str) -> str:
    if os.path.sep in command:
        return command

    sibling = Path(sys.executable).resolve().parent / command
    if sibling.is_file() and os.access(sibling, os.X_OK):
        return str(sibling)

    resolved = shutil.which(command)
    if resolved:
        return resolved

    raise SystemExit(
        f"Missing {command!r}. Run this script with the pySCENIC conda Python "
        "or pass --pyscenic-command /path/to/pyscenic."
    )


def run_command(label: str, command: list[str], log_path: Path | None = None) -> None:
    print(label, flush=True)
    print("[CMD] " + " ".join(command), flush=True)
    if log_path is None:
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as exc:
            raise SystemExit(f"{label} failed with exit status {exc.returncode}.") from exc
        return

    log_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[LOG] {log_path}", flush=True)
    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write(label + "\n")
        log_handle.write("[CMD] " + " ".join(command) + "\n\n")
        log_handle.flush()
        try:
            subprocess.run(command, stdout=log_handle, stderr=subprocess.STDOUT, check=True)
        except subprocess.CalledProcessError as exc:
            raise SystemExit(
                f"{label} failed with exit status {exc.returncode}. See log: {log_path}"
            ) from exc


def main() -> int:
    args = parse_args()
    configure_subprocess_threads()

    # Validate dependencies and biological identifier compatibility before
    # starting expensive network inference.
    validate_pyscenic_python_environment()

    check_file(args.expr_csv, "expression matrix")
    check_file(args.tf_list, "TF list")
    check_file(args.motif_annotations, "motif annotations")
    for ranking_db in args.ranking_db:
        check_file(ranking_db, "ranking database")
    validate_tf_overlap(args.expr_csv, args.tf_list)

    pyscenic = resolve_pyscenic_command(args.pyscenic_command)
    log_dir = Path(args.log_dir) if args.log_dir else None
    grn_num_workers = args.grn_num_workers or default_grn_workers(args.num_workers)
    parse_positive_int(grn_num_workers, "--grn-num-workers")
    expr_loom = args.expr_loom or str(
        Path(args.auc_loom).with_name("expression_matrix.loom")
    )

    # Stage 1: infer TF-target relationships. GRN worker count defaults to one
    # in the shell wrapper because local Dask workers were unstable at scale.
    grn_command = [
        pyscenic,
        "grn",
        args.expr_csv,
        args.tf_list,
        "--output",
        args.adj_tsv,
        "--method",
        args.grn_method,
        "--num_workers",
        str(grn_num_workers),
    ]
    if args.seed:
        grn_command.extend(["--seed", str(args.seed)])

    run_command(
        "Running pySCENIC GRN inference",
        grn_command,
        log_dir / "pyscenic_grn.log" if log_dir else None,
    )
    validate_grn_output(
        args.adj_tsv,
        args.expr_csv,
        args.tf_list,
        log_dir / "pyscenic_grn.log" if log_dir else None,
    )

    # Stage 2: retain motif-supported TF-target modules using the cisTarget
    # ranking databases and motif annotation table.
    run_command(
        "Running pySCENIC cisTarget context pruning",
        [
            pyscenic,
            "ctx",
            args.adj_tsv,
            *args.ranking_db,
            "--annotations_fname",
            args.motif_annotations,
            "--expression_mtx_fname",
            args.expr_csv,
            "--mode",
            args.ctx_mode,
            "--output",
            args.regulons_csv,
            "--num_workers",
            str(args.num_workers),
        ],
        log_dir / "pyscenic_ctx.log" if log_dir else None,
    )

    # Stage 3: build or validate the expression loom needed by AUCell. GRN and
    # ctx intentionally continue using CSV because those commands accept it.
    expr_loom_path = Path(expr_loom)
    expr_csv_path = Path(args.expr_csv)
    if (
        not expr_loom_path.is_file()
        or expr_loom_path.stat().st_size == 0
        or expr_loom_path.stat().st_mtime < expr_csv_path.stat().st_mtime
    ):
        convert_expression_csv_to_loom(
            args.expr_csv,
            expr_loom,
            args.loom_chunk_cells,
        )
    else:
        print(f"[INFO] Reusing expression loom: {expr_loom}", flush=True)
        try:
            validate_expression_loom(expr_loom)
        except SystemExit:
            print(
                f"[WARN] Existing expression loom is invalid; rebuilding: {expr_loom}",
                flush=True,
            )
            convert_expression_csv_to_loom(
                args.expr_csv,
                expr_loom,
                args.loom_chunk_cells,
            )

    # Stage 4: score regulon activity per cell. Remove stale output first so a
    # failed command cannot leave an old auc_mtx.loom that appears successful.
    auc_log = log_dir / "pyscenic_aucell.log" if log_dir else None
    auc_path = Path(args.auc_loom)
    if auc_path.exists():
        auc_path.unlink()

    run_command(
        "Running pySCENIC AUCell scoring",
        [
            pyscenic,
            "aucell",
            expr_loom,
            args.regulons_csv,
            "--output",
            args.auc_loom,
            "--num_workers",
            str(args.num_workers),
        ],
        auc_log,
    )
    # A zero exit status is insufficient: validate both HDF5 structure and the
    # pySCENIC-specific CellID/RegulonsAUC datasets.
    validate_auc_loom(args.auc_loom, auc_log)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
