#!/usr/bin/env python3
"""Semantic integrity audit for the frozen T1/T2/T_R confirmation."""
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


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path: Path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def num(row, key):
    return float(row[key])


def integer(row, key):
    return int(row[key])


def summarize(group):
    err = np.asarray([num(r, "q_error_riu") for r in group])
    sq = err**2
    n = len(err)
    return {
        "n": n,
        "bias_riu": float(np.mean(err)),
        "bias_mcse_riu": float(np.std(err, ddof=1) / math.sqrt(n)),
        "variance_riu2": float(np.var(err, ddof=1)),
        "mse_riu2": float(np.mean(sq)),
        "mse_mcse_riu2": float(np.std(sq, ddof=1) / math.sqrt(n)),
        "rmse_riu": float(np.sqrt(np.mean(sq))),
        "mae_riu": float(np.mean(np.abs(err))),
        "p95_absolute_error_riu": float(np.quantile(np.abs(err), 0.95)),
        "optimizer_failure_rate": float(np.mean([integer(r, "optimizer_failure") for r in group])),
        "target_boundary_rate": float(np.mean([integer(r, "target_boundary") for r in group])),
        "multi_solution_rate": float(np.mean([integer(r, "multi_solution") for r in group])),
        "median_final_fisher_step_inf_norm_unit": float(np.median([num(r, "final_fisher_step_inf_norm_unit") for r in group])),
        "max_final_fisher_step_inf_norm_unit": float(np.max([num(r, "final_fisher_step_inf_norm_unit") for r in group])),
    }


def main(config_path):
    cp = Path(config_path).resolve()
    cfg = json.loads(cp.read_text(encoding="utf-8"))
    out = ROOT / "outputs" / "integrity_audit" / cfg["run_id"]
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True)

    paths = {k: ROOT / cfg[k] for k in ("input_replicates", "input_summary", "frozen_config", "generator_script")}
    for name, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"missing {name}: {path}")
    hashes = {name: digest(path) for name, path in paths.items()}
    rows = read_csv(paths["input_replicates"])
    frozen = json.loads(paths["frozen_config"].read_text(encoding="utf-8"))
    reported = json.loads(paths["input_summary"].read_text(encoding="utf-8"))

    required = {
        "task_id", "scenario_id", "scenario", "replicate_id", "strategy", "fit_stage",
        "config_hash", "noise_stream_id", "true_mismatch_nm", "true_target_angle_offset_deg",
        "source_reference_fraction", "target_reference_fraction", "target_sample_fraction",
        "truth_q", "estimated_q", "q_error_riu", "objective_half_deviance",
        "optimizer_failure", "target_boundary", "multi_solution",
        "final_fisher_step_inf_norm_unit",
    }
    missing_columns = sorted(required - set(rows[0])) if rows else sorted(required)
    unique_keys = [(integer(r, "scenario_id"), integer(r, "replicate_id"), r["strategy"]) for r in rows]
    duplicate_keys = [list(k) + [n] for k, n in Counter(unique_keys).items() if n != 1]

    expected = set()
    expected_task_ids = set()
    condition_by_id = {int(c["condition_id"]): c for c in frozen["conditions"]}
    for c in frozen["conditions"]:
        cid = int(c["condition_id"])
        for rep in range(int(frozen["replicates_per_condition"])):
            for strategy in ("T1", "T2", "T_R"):
                expected.add((cid, rep, strategy))
                expected_task_ids.add(f"c{cid}:r{rep}:{strategy}")
    observed = set(unique_keys)
    missing_tasks = sorted([list(x) for x in expected - observed])
    orphan_rows = sorted([list(x) for x in observed - expected])
    task_id_mismatch = sorted(r["task_id"] for r in rows if r["task_id"] not in expected_task_ids)
    wrong_stage = sorted(r["task_id"] for r in rows if r["fit_stage"] != "final_refined")
    wrong_row_config_hash = sorted(r["task_id"] for r in rows if r["config_hash"] != hashes["frozen_config"])

    allocation_mismatches = []
    stream_violations = []
    grouped_streams = defaultdict(dict)
    for r in rows:
        cid = integer(r, "scenario_id")
        rep = integer(r, "replicate_id")
        strategy = r["strategy"]
        c = condition_by_id.get(cid)
        if c is None or r["scenario"] != c["scenario"] or num(r, "true_mismatch_nm") != float(c["mismatch_nm"]) or num(r, "true_target_angle_offset_deg") != float(c["target_angle_offset_deg"]):
            allocation_mismatches.append([r["task_id"], "condition_metadata"])
            continue
        mode = "angle_aware" if bool(c["angle_aware"]) else "stable_angle"
        wanted = frozen["allocations"][mode][strategy]
        got = [num(r, "source_reference_fraction"), num(r, "target_reference_fraction"), num(r, "target_sample_fraction")]
        if any(abs(a - b) > 1e-12 for a, b in zip(got, wanted)):
            allocation_mismatches.append([r["task_id"], "allocation"])
        grouped_streams[(cid, rep)][strategy] = r["noise_stream_id"]
    for key, streams in grouped_streams.items():
        if set(streams) != {"T1", "T2", "T_R"} or len(set(streams.values())) != 3:
            stream_violations.append([key[0], key[1], json.dumps(streams, sort_keys=True)])

    groups = defaultdict(list)
    for r in rows:
        groups[(integer(r, "scenario_id"), r["strategy"])].append(r)
    verified = []
    count_mismatches = []
    for key in sorted(groups):
        item = {"scenario_id": key[0], "scenario": groups[key][0]["scenario"], "strategy": key[1], **summarize(groups[key])}
        verified.append(item)
        if item["n"] != int(frozen["replicates_per_condition"]):
            count_mismatches.append([key[0], key[1], item["n"]])

    by = {(r["scenario_id"], r["strategy"]): r for r in verified}
    differences = []
    for cid in sorted(condition_by_id):
        for a, b in (("T1", "T_R"), ("T2", "T_R"), ("T1", "T2")):
            ra, rb = by[(cid, a)], by[(cid, b)]
            mean = ra["mse_riu2"] - rb["mse_riu2"]
            se = math.sqrt(ra["mse_mcse_riu2"]**2 + rb["mse_mcse_riu2"]**2)
            differences.append({
                "scenario_id": cid, "scenario": ra["scenario"], "strategy_a": a, "strategy_b": b,
                "comparison_design": "independent protocol-specific Poisson streams; not paired",
                "mean_mse_difference_a_minus_b_riu2": mean, "mcse_riu2": se,
                "mc_interval_low_riu2": mean - 1.96 * se,
                "mc_interval_high_riu2": mean + 1.96 * se,
            })

    optimizer_failures = sum(integer(r, "optimizer_failure") for r in rows)
    target_boundary_hits = sum(integer(r, "target_boundary") for r in rows)
    multi_solution_flags = sum(integer(r, "multi_solution") for r in rows)
    checks = {
        "record_count": len(rows) == int(cfg["expected_records"]),
        "unique_expected_tasks": not duplicate_keys and not missing_tasks and not orphan_rows,
        "task_ids": not task_id_mismatch,
        "required_columns": not missing_columns,
        "final_stage_only": not wrong_stage,
        "frozen_config_hash": hashes["frozen_config"] == cfg["expected_frozen_config_sha256"],
        "generator_hash": hashes["generator_script"] == cfg["expected_generator_script_sha256"],
        "row_config_hashes": not wrong_row_config_hash,
        "condition_and_allocation_metadata": not allocation_mismatches,
        "independent_method_streams": not stream_violations,
        "cell_counts": len(verified) == 12 and not count_mismatches,
        "reported_status": reported.get("status") == "FROZEN_TARGET_RECALIBRATION_CONFIRMATION_COMPLETE",
        "reported_counts": reported.get("task_count") == int(cfg["expected_records"]) and reported.get("record_count") == int(cfg["expected_records"]),
    }
    report = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "run_id": cfg["run_id"], "audit_seed": cfg["seed"], "checks": checks,
        "hashes": hashes, "reported_status": reported.get("status"),
        "record_count": len(rows), "expected_records": cfg["expected_records"],
        "cell_count": len(verified), "expected_cells": 12,
        "missing_columns": missing_columns, "duplicate_unique_keys": duplicate_keys,
        "missing_tasks": missing_tasks, "orphan_rows": orphan_rows,
        "task_id_mismatches": task_id_mismatch, "wrong_fit_stage": wrong_stage,
        "wrong_row_config_hash": wrong_row_config_hash,
        "allocation_or_condition_mismatches": allocation_mismatches,
        "method_stream_independence_violations": stream_violations,
        "count_mismatches": count_mismatches,
        "optimizer_failures": optimizer_failures,
        "target_boundary_hits": target_boundary_hits,
        "multi_solution_flags": multi_solution_flags,
        "comparison_note": "Cross-method risk differences use independent Monte Carlo standard errors because each method used a distinct protocol-specific Poisson stream.",
    }
    write_csv(out / "condition_summary_verified.csv", verified)
    write_csv(out / "independent_risk_differences_verified.csv", differences)
    (out / "integrity_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)
    if report["status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_target_recalibration_audit.py CONFIG")
    main(sys.argv[1])
