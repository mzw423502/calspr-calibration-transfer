#!/usr/bin/env python3
"""Semantic completeness and paired-risk audit for the frozen confirmation."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def f(row, key):
    return float(row[key])


def i(row, key):
    return int(row[key])


def summarize(group):
    err = np.asarray([f(r, "q_error_riu") for r in group])
    sq = err**2
    n = len(err)
    return {
        "n": n,
        "bias_riu": float(np.mean(err)),
        "bias_mcse_riu": float(np.std(err, ddof=1)/math.sqrt(n)),
        "variance_riu2": float(np.var(err, ddof=1)),
        "mse_riu2": float(np.mean(sq)),
        "mse_mcse_riu2": float(np.std(sq, ddof=1)/math.sqrt(n)),
        "rmse_riu": float(np.sqrt(np.mean(sq))),
        "mae_riu": float(np.mean(np.abs(err))),
        "p95_absolute_error_riu": float(np.quantile(np.abs(err), 0.95)),
        "max_absolute_error_riu": float(np.max(np.abs(err))),
        "optimizer_failure_rate": float(np.mean([i(r, "optimizer_failure") for r in group])),
        "target_boundary_rate": float(np.mean([i(r, "target_boundary") for r in group])),
        "multi_solution_rate": float(np.mean([i(r, "multi_solution") for r in group])),
        "median_final_fisher_step_inf_norm_unit": float(np.median([f(r, "final_fisher_step_inf_norm_unit") for r in group])),
        "max_final_fisher_step_inf_norm_unit": float(np.max([f(r, "final_fisher_step_inf_norm_unit") for r in group])),
    }


def paired_summary(a, b, label_a, label_b, scenario, delta, tau):
    map_a = {(i(r, "replicate"), r["noise_stream_id"]): r for r in a}
    map_b = {(i(r, "replicate"), r["noise_stream_id"]): r for r in b}
    keys = sorted(set(map_a) & set(map_b))
    d = np.asarray([f(map_a[k], "q_error_riu")**2-f(map_b[k], "q_error_riu")**2 for k in keys])
    n = len(d)
    mean = float(np.mean(d))
    se = float(np.std(d, ddof=1)/math.sqrt(n))
    return {
        "scenario": scenario, "true_mismatch_nm": delta,
        "true_target_angle_offset_deg": tau,
        "strategy_a": label_a, "strategy_b": label_b,
        "paired_key": "replicate+noise_stream_id", "n_pairs": n,
        "mean_mse_difference_a_minus_b_riu2": mean,
        "mcse_riu2": se, "mc_interval_low_riu2": mean-1.96*se,
        "mc_interval_high_riu2": mean+1.96*se,
        "fraction_difference_positive": float(np.mean(d > 0)),
    }


def main(config_path):
    cfg_path = Path(config_path).resolve()
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    out = ROOT / "outputs" / "integrity_audit" / cfg["run_id"]
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True)

    paths = {k: ROOT/cfg[k] for k in ("input_replicates", "input_summary", "input_condition_summary", "frozen_config", "generator_script", "boundary_exploration_summary")}
    for name, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"missing {name}: {path}")
    hashes = {name: digest(path) for name, path in paths.items()}
    rows = read_csv(paths["input_replicates"])
    frozen = json.loads(paths["frozen_config"].read_text(encoding="utf-8"))
    run_summary = json.loads(paths["input_summary"].read_text(encoding="utf-8"))
    boundary = json.loads(paths["boundary_exploration_summary"].read_text(encoding="utf-8"))

    required = {"scenario", "condition_id", "strategy", "replicate", "noise_stream_id",
                "true_mismatch_nm", "true_target_angle_offset_deg", "truth_q", "estimated_q",
                "q_error_riu", "objective_half_deviance", "optimizer_failure", "target_boundary",
                "multi_solution", "final_fisher_step_inf_norm_unit"}
    missing_columns = sorted(required-set(rows[0])) if rows else sorted(required)
    keys = [(r["scenario"], i(r, "condition_id"), r["strategy"], i(r, "replicate")) for r in rows]
    duplicate_keys = [list(k)+[v] for k, v in Counter(keys).items() if v != 1]

    expected = []
    cid = 0
    for delta in frozen["mismatch_nm"]:
        expected.append((cid, "thickness_mismatch", float(delta), 0.0,
                         {"pair": ["T0", "T1", "O"], "t2": ["T2"]})); cid += 1
    for tau in frozen["target_angle_offset_deg"]:
        expected.append((cid, "target_angle_offset", 0.0, float(tau),
                         {"pair": ["T1", "T1_angle"], "t2": ["T2", "T2_angle"]})); cid += 1

    by_condition_strategy = defaultdict(list)
    for row in rows:
        by_condition_strategy[(i(row, "condition_id"), row["strategy"])].append(row)
    manifest = []
    missing_tasks = []
    orphan_rows = []
    expected_row_keys = set()
    for condition_id, scenario, delta, tau, task_types in expected:
        for rep in range(frozen["replicates_per_condition"]):
            for task_type, strategies in task_types.items():
                observed = []
                streams = set()
                for strategy in strategies:
                    key = (scenario, condition_id, strategy, rep)
                    expected_row_keys.add(key)
                    matches = [r for r in by_condition_strategy[(condition_id, strategy)] if i(r, "replicate") == rep]
                    if len(matches) == 1:
                        observed.append(strategy); streams.add(matches[0]["noise_stream_id"])
                complete = observed == strategies and len(streams) == 1
                task_id = f"c{condition_id}:r{rep}:{task_type}"
                manifest.append({"task_id":task_id, "scenario_id":condition_id, "scenario":scenario,
                                 "replicate_id":rep, "task_type":task_type,
                                 "expected_strategies":"|".join(strategies),
                                 "observed_strategies":"|".join(observed),
                                 "shared_noise_stream_within_task":int(len(streams) == 1),
                                 "complete":int(complete), "config_hash":hashes["frozen_config"]})
                if not complete:
                    missing_tasks.append(task_id)
    for row in rows:
        key = (row["scenario"], i(row, "condition_id"), row["strategy"], i(row, "replicate"))
        if key not in expected_row_keys:
            orphan_rows.append(list(key))

    summaries = []
    groups = defaultdict(list)
    for row in rows:
        groups[(row["scenario"], f(row, "true_mismatch_nm"), f(row, "true_target_angle_offset_deg"), row["strategy"])].append(row)
    for key in sorted(groups, key=lambda z: tuple(str(x) for x in z)):
        summaries.append({"scenario":key[0], "true_mismatch_nm":key[1],
                          "true_target_angle_offset_deg":key[2], "strategy":key[3],
                          **summarize(groups[key])})

    paired = []
    for condition_id, scenario, delta, tau, _ in expected:
        if scenario == "thickness_mismatch":
            paired.append(paired_summary(by_condition_strategy[(condition_id, "T0")], by_condition_strategy[(condition_id, "T1")], "T0", "T1", scenario, delta, tau))
        else:
            paired.append(paired_summary(by_condition_strategy[(condition_id, "T1")], by_condition_strategy[(condition_id, "T1_angle")], "T1", "T1_angle", scenario, delta, tau))
            paired.append(paired_summary(by_condition_strategy[(condition_id, "T2")], by_condition_strategy[(condition_id, "T2_angle")], "T2", "T2_angle", scenario, delta, tau))

    thickness_pairs = [r for r in paired if r["scenario"] == "thickness_mismatch"]
    pooled = []
    for magnitude in sorted({abs(float(r["true_mismatch_nm"])) for r in thickness_pairs}):
        source = [r for r in thickness_pairs if abs(float(r["true_mismatch_nm"])) == magnitude]
        raw_d = []
        for s in source:
            condition_id = next(x[0] for x in expected if x[1] == "thickness_mismatch" and x[2] == s["true_mismatch_nm"])
            a = by_condition_strategy[(condition_id, "T0")]
            b = by_condition_strategy[(condition_id, "T1")]
            ma = {(i(r,"replicate"),r["noise_stream_id"]):r for r in a}
            mb = {(i(r,"replicate"),r["noise_stream_id"]):r for r in b}
            raw_d.extend(f(ma[k],"q_error_riu")**2-f(mb[k],"q_error_riu")**2 for k in set(ma)&set(mb))
        arr = np.asarray(raw_d); mean = float(np.mean(arr)); se = float(np.std(arr,ddof=1)/math.sqrt(len(arr)))
        pooled.append({"absolute_mismatch_nm":magnitude, "n_pairs":len(arr),
                       "mean_mse_difference_T0_minus_T1_riu2":mean, "mcse_riu2":se,
                       "mc_interval_low_riu2":mean-1.96*se, "mc_interval_high_riu2":mean+1.96*se})
    negative = [r for r in pooled if r["mean_mse_difference_T0_minus_T1_riu2"] < 0]
    positive = [r for r in pooled if r["mean_mse_difference_T0_minus_T1_riu2"] > 0]
    empirical_bracket = {
        "largest_tested_magnitude_with_negative_mean_D_nm": max((r["absolute_mismatch_nm"] for r in negative), default=None),
        "smallest_tested_magnitude_with_positive_mean_D_nm": min((r["absolute_mismatch_nm"] for r in positive), default=None),
        "predicted_local_boundary_nm": boundary["anchor_T0_vs_T1_boundary_abs_delta_nm"],
        "interpretation": "A sign bracket from the tested grid, not a fitted threshold or a device-population interval. Intervals crossing zero indicate limited Monte Carlo resolution near the boundary."
    }

    write_csv(out/"task_manifest_verified.csv", manifest)
    write_csv(out/"condition_summary_verified.csv", summaries)
    write_csv(out/"paired_risk_differences.csv", paired)
    write_csv(out/"symmetric_magnitude_risk_differences.csv", pooled)
    report = {
        "status": "PASS" if (not missing_columns and not duplicate_keys and not missing_tasks and not orphan_rows and len(rows) == cfg["expected_records"] and len(manifest) == cfg["expected_tasks"] and hashes["frozen_config"] == cfg["expected_frozen_config_sha256"] and hashes["generator_script"] == cfg["expected_generator_script_sha256"] and run_summary.get("status") == "FROZEN_INTERNAL_CONFIRMATION_COMPLETE") else "FAIL",
        "run_id":cfg["run_id"], "audit_seed":cfg["seed"], "hashes":hashes,
        "generator_reported_status":run_summary.get("status"),
        "record_count":len(rows), "expected_records":cfg["expected_records"],
        "task_count":len(manifest), "expected_tasks":cfg["expected_tasks"],
        "missing_columns":missing_columns, "duplicate_unique_keys":duplicate_keys,
        "missing_or_incomplete_tasks":missing_tasks, "orphan_rows":orphan_rows,
        "schema_note":"The generator did not store task_id, fit_stage, or config_hash per row. This audit reconstructs one final-refined row per strategy and emits a task manifest with the frozen config hash; it does not relabel initial and refined fits as independent replicates.",
        "interval_method":cfg["interval_method"], "empirical_crossing_assessment":empirical_bracket,
        "optimizer_failures":sum(i(r,"optimizer_failure") for r in rows),
        "target_boundary_hits":sum(i(r,"target_boundary") for r in rows),
        "multi_solution_flags":sum(i(r,"multi_solution") for r in rows),
    }
    (out/"integrity_report.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)
    if report["status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_boundary_confirmation_audit.py CONFIG")
    main(sys.argv[1])
