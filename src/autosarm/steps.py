"""Stateless CLI step for Argo Workflows (ADR 0012 P1/P2).

    python -m autosarm.steps run --work-dir /work    # CPU only

⚠️ **This has never run.** There is no `autosarm` image on the node and no Deployment in
the cluster, so `sarm_analysis` tasks have never been executed by anything. It is written
to the same contract as the other steps and should work once the image is built, but do
not mistake it for verified code.

Replaces `tasks/processor.py`. Of the four workers, autoSARM's polling was the only
*correct* one (`FOR UPDATE SKIP LOCKED` inside a real transaction) — and it is deleted
anyway, because Argo owns the claim now and there is nothing left to lock.

The science layer is already CLI-shaped (`autosarm sarm` / `autosarm tree` build an
argparse.Namespace and call `run_sarm` / `run_tree`), so this step does the same thing the
worker did: synthesise the Namespace from params.json and call in.

## The subtype dispatch

`sarm_analysis` is really two different jobs sharing one task_type, told apart only by
`sarm_task_params.task_subtype` ('sarm' | 'tree'). We dispatch on it here rather than
building two WorkflowTemplates, because the backend gives us no way to tell them apart at
submit time without reading that same column.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from argparse import Namespace
from pathlib import Path
from typing import Any, Dict

logger = logging.getLogger("sarm-step")


def _run_matrix(work_dir: Path, params: Dict[str, Any]) -> None:
    from autosarm.cli.create_sarm import run_sarm

    csv_name = params.get("csv_filename", "compounds.csv")
    csv_path = work_dir / "input" / csv_name
    if not csv_path.exists():
        raise FileNotFoundError(f"input CSV {csv_name} not found in the job's input/")

    value_columns = params.get("value_columns") or "[]"
    if isinstance(value_columns, str):
        value_columns = json.loads(value_columns)

    run_sarm(
        Namespace(
            csvFile=str(csv_path),
            type=params.get("analysis_type", "smiles"),
            column=value_columns,
            log=bool(params.get("log_transform", False)),
            minimumSite1=int(params.get("minimum_site1", 3)),
            minimumSite2=int(params.get("minimum_site2", 3)),
            n_jobs=int(params.get("n_jobs", 8)),
            save_folder=str(work_dir / "output" / "SAR_Results"),
            csv2excel=bool(params.get("csv2excel", False)),
        )
    )


def _run_tree(work_dir: Path, params: Dict[str, Any]) -> None:
    from autosarm.cli.create_tree import run_tree

    core = params.get("fragment_core")
    if not core:
        raise ValueError("tree subtype requires `fragment_core`")

    # tree.py expects Combine_Table_info.csv etc. directly in work_folder, and the matrix
    # stage puts them inside a SAR_Results/ subdirectory — so point it at that subdirectory.
    # The backend has already copied the upstream SAR_Results into this job's own output/
    # prefix (routers/sarm.py:380-384), and `fetch` does not touch output/, so we recreate
    # that expectation locally instead.
    work_folder = work_dir / "output" / "SAR_Results"
    if not (work_folder / "Combine_Table_info.csv").exists():
        raise FileNotFoundError(
            "no Combine_Table_info.csv — the upstream SAR matrix results were not copied "
            "into this job (the backend does that at submit time from source_task_id)"
        )

    run_tree(
        Namespace(
            fragment_core=core,
            rootTitle=params.get("root_title", "root"),
            workFolder=str(work_folder),
            inputFile=params.get("input_file", "input.csv"),
            maxLevel=int(params.get("max_level", 5)),
            treeContent=str(params.get("tree_content", '["double-cut"]')),
            highlightDict=str(params.get("highlight_dict", "[]")),
        )
    )


def stage_run(work_dir: Path, params: Dict[str, Any]) -> None:
    subtype = params.get("task_subtype", "sarm")
    (work_dir / "output").mkdir(parents=True, exist_ok=True)
    logger.info("sarm subtype=%s", subtype)
    if subtype == "sarm":
        _run_matrix(work_dir, params)
    elif subtype == "tree":
        _run_tree(work_dir, params)
    else:
        raise ValueError(f"unknown task_subtype {subtype!r} (expected 'sarm' or 'tree')")


STAGES = {"run": stage_run}


def main() -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stdout,
    )
    p = argparse.ArgumentParser(prog="autosarm.steps")
    p.add_argument("stage", choices=sorted(STAGES))
    p.add_argument("--work-dir", default="/work")
    args = p.parse_args()

    work_dir = Path(args.work_dir)
    params = json.loads((work_dir / "params.json").read_text())

    try:
        STAGES[args.stage](work_dir, params)
    except Exception as e:
        print(f"FATAL [{args.stage}] {type(e).__name__}: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
