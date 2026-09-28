#!/usr/bin/env python3
"""Run the predeclared CalSPR V2.4 precision extension for replicates 100-199."""
from __future__ import annotations

import csv
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import run_v24_confirmation as base


def main(config_path):
    started = time.time()
    cfg_path = Path(config_path).resolve()
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    source = ROOT / cfg["selection_source"]
    trigger = ROOT / cfg["extension_trigger_source"]
    if base.digest(source) != cfg["selection_source_sha256"]:
        raise RuntimeError("threshold prediction source hash mismatch")
    if base.digest(trigger) != cfg["extension_trigger_source_sha256"]:
        raise RuntimeError("extension trigger source hash mismatch")
    trigger_summary = json.loads(trigger.read_text(encoding="utf-8"))
    if not trigger_summary.get("extension_rule_triggered"):
        raise RuntimeError("frozen source does not trigger the extension rule")
    start = int(cfg["replicate_start"])
    stop = start + int(cfg["replicates_per_condition"])
    if start != 100 or stop != 200:
        raise RuntimeError("this frozen extension must cover replicates 100-199")
    out = ROOT / "outputs" / "v24" / cfg["run_id"]
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True)
    tasks = [(cfg, anchor, float(rho), rep) for anchor in cfg["anchors"]
             for rho in cfg["relative_magnitudes"] for rep in range(start, stop)]
    partial = out / "V24_REPLICATES.partial.csv"
    rows = []
    with open(partial, "w", newline="", encoding="utf-8") as handle, ProcessPoolExecutor(max_workers=cfg["worker_processes"]) as pool:
        writer = None
        futures = [pool.submit(base.paired_task, task) for task in tasks]
        for done, future in enumerate(as_completed(futures), start=1):
            pair = future.result()
            if writer is None:
                writer = csv.DictWriter(handle, fieldnames=list(pair[0]))
                writer.writeheader()
            writer.writerows(pair); handle.flush(); rows.extend(pair)
            if done % 100 == 0:
                print(f"completed extension paired tasks {done}/{len(tasks)}", flush=True)
    rows.sort(key=lambda r: (r["anchor_id"], float(r["rho_delta_over_pred"]), r["replicate"], r["strategy"]))
    final_replicates = out / "V24_REPLICATES.csv"
    base.write_csv(final_replicates, rows)
    partial.unlink()
    summaries, extend = base.summarize(rows, cfg)
    base.write_csv(out / "V24_CONFIRMATION_SUMMARY.csv", summaries)
    expected = len(tasks) * 2
    summary = {
        "status": "V24_CONFIRMATION_EXTENSION_COMPLETE",
        "run_id": cfg["run_id"], "seed": cfg["seed"],
        "replicate_start": start, "replicate_stop_exclusive": stop,
        "replicates_per_condition": cfg["replicates_per_condition"],
        "worker_processes": cfg["worker_processes"], "task_count": len(tasks),
        "record_count": len(rows), "expected_record_count": expected,
        "record_count_pass": len(rows) == expected,
        "optimizer_failure_count": int(sum(r["optimizer_failure"] for r in rows)),
        "target_boundary_count": int(sum(r["target_boundary"] for r in rows)),
        "multi_solution_count": int(sum(r["multi_solution"] for r in rows)),
        "extension_rule_still_triggered_on_increment_only": bool(extend),
        "elapsed_seconds": time.time() - started,
        "config_sha256": base.digest(cfg_path),
        "selection_source_sha256": base.digest(source),
        "extension_trigger_source_sha256": base.digest(trigger),
        "replicates_sha256": base.digest(final_replicates),
        "interpretation_limit": cfg["interpretation_limit"],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_v24_confirmation_extension.py CONFIG")
    main(sys.argv[1])
