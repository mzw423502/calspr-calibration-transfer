#!/usr/bin/env python3
"""Merge the frozen n=100 run and its predeclared 100-199 extension."""
from __future__ import annotations

import csv
import json
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import run_v24_confirmation as base

BASE = ROOT / "outputs/v24/v24_confirmation_20260914_c1_fast_validated"
EXT = ROOT / "outputs/v24/v24_confirmation_20260914_c2_extension_rep100_199"
OUT = ROOT / "outputs/v24/v24_confirmation_20260914_c3_combined_n200"
CONFIG = ROOT / "configs/v24_confirmation_extension_v1.json"

FLOATS = {"thickness_nm","total_incident_photons","rho_delta_over_pred","delta_pred_nm","true_delta_nm","true_q","error_riu","squared_error_riu2","estimated_q","objective_half_deviance","final_score_scaled_inf","final_fisher_step_inf_norm_unit"}
INTS = {"replicate","optimizer_failure","target_boundary","multi_solution","fit_iterations","fit_accepted_steps"}


def load(path):
    with open(path, newline="", encoding="utf-8") as f:
        data = list(csv.DictReader(f))
    for row in data:
        for key in FLOATS: row[key] = float(row[key])
        for key in INTS: row[key] = int(row[key])
    return data


def main():
    started=time.time(); cfg=json.loads(CONFIG.read_text())
    if OUT.exists(): raise FileExistsError(f"refusing to overwrite {OUT}")
    required={
      BASE/"V24_REPLICATES.csv":"3ff79b74c03c6c3856afed78bd8075df4c78bebc19a0884e25e663595cb11417",
      EXT/"V24_REPLICATES.csv":"ed9a1e7e37d60dc4a2ef6697f1dd573d1f3f50f4ae9cf6c9ce84d8f8c9d5fa61",
    }
    for path, expected in required.items():
        got=base.digest(path)
        if got != expected: raise RuntimeError(f"source hash mismatch: {path}: {got}")
    rows=load(BASE/"V24_REPLICATES.csv")+load(EXT/"V24_REPLICATES.csv")
    keys=[(r["task_id"],r["strategy"]) for r in rows]
    if len(rows)!=25200 or len(set(keys))!=25200: raise RuntimeError("record count or uniqueness failure")
    reps={int(r["replicate"]) for r in rows}
    if reps != set(range(200)): raise RuntimeError("replicate coverage is not exactly 0-199")
    counts=Counter((r["anchor_id"],r["rho_delta_over_pred"],r["strategy"]) for r in rows)
    if len(counts)!=126 or set(counts.values())!={200}: raise RuntimeError("condition counts are not 126 groups of n=200")
    rows.sort(key=lambda r:(r["anchor_id"],r["rho_delta_over_pred"],r["replicate"],r["strategy"]))
    OUT.mkdir(parents=True)
    final=OUT/"V24_REPLICATES.csv"; base.write_csv(final,rows)
    summaries,trigger=base.summarize(rows,cfg); base.write_csv(OUT/"V24_CONFIRMATION_SUMMARY.csv",summaries)
    diagnostics={k:int(sum(r[k] for r in rows)) for k in ("optimizer_failure","target_boundary","multi_solution")}
    if any(diagnostics.values()): raise RuntimeError(f"nonzero diagnostics: {diagnostics}")
    summary={
      "status":"V24_CONFIRMATION_COMBINED_N200_COMPLETE",
      "run_id":"v24_confirmation_20260914_c3_combined_n200",
      "seed":cfg["seed"],"replicates_per_condition":200,"replicate_range":[0,199],
      "task_count":12600,"record_count":len(rows),"expected_record_count":25200,
      "condition_strategy_group_count":len(counts),"all_group_n":sorted(set(counts.values())),
      "optimizer_failure_count":diagnostics["optimizer_failure"],
      "target_boundary_count":diagnostics["target_boundary"],
      "multi_solution_count":diagnostics["multi_solution"],
      "precision_rule_triggered_after_maximum_n200":bool(trigger),
      "precision_extension_limit_reached":True,
      "further_extension_prohibited_by_frozen_rule":True,
      "elapsed_seconds":time.time()-started,
      "base_replicates_sha256":required[BASE/"V24_REPLICATES.csv"],
      "extension_replicates_sha256":required[EXT/"V24_REPLICATES.csv"],
      "combined_replicates_sha256":base.digest(final),
      "combined_summary_sha256":base.digest(OUT/"V24_CONFIRMATION_SUMMARY.csv"),
      "interpretation_limit":cfg["interpretation_limit"],
    }
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__": main()
