#!/usr/bin/env python3
"""Frozen noisy comparison of T1, T2, and target-only recalibration T_R."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"scripts"))
from calspr.inversion import bounded_pattern_fit, fisher_scoring_refine, poisson_deviance_half
from run_boundary_confirmation import final_diagnostics, fit_one
from run_boundary_scan import model_vector, truth_vector
from run_target_recalibration_scan import tr_bounds, tr_model, tr_truth


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fit_tr(y, cfg, angles, gold, target_reference_fraction):
    bounds = tr_bounds(cfg)
    fn = lambda p: tr_model(p, cfg, angles, gold, target_reference_fraction)
    sols = bounded_pattern_fit(lambda p: poisson_deviance_half(y, fn(p)), bounds, cfg["fit"])
    best = sols[0]
    refined = fisher_scoring_refine(y, fn, best["params"], bounds,
                                    max_iter=cfg["fit"]["refine_max_iterations"],
                                    derivative_step=cfg["fit"]["refine_unit_derivative_step"])
    target = float(refined["params"][0])
    lo, hi = bounds[0]
    bf = cfg["fit"]["target_boundary_fraction"]
    multi = any(s["objective"]-best["objective"] <= cfg["fit"]["multi_solution_deviance_delta"]
                and abs(float(s["params"][0])-float(best["params"][0])) >= cfg["fit"]["multi_solution_target_separation_riu"]
                for s in sols[1:])
    score_inf, step_inf = final_diagnostics(y, fn, refined["params"], bounds,
                                            cfg["fit"]["refine_unit_derivative_step"])
    return {
        "estimated_q":target, "objective_half_deviance":float(refined["objective"]),
        "optimizer_failure":int(not np.isfinite(refined["objective"])),
        "target_boundary":int(target-lo <= bf*(hi-lo) or hi-target <= bf*(hi-lo)),
        "multi_solution":int(multi), "refine_iterations":refined["iterations"],
        "refine_accepted_steps":refined["accepted_steps"],
        "final_score_inf_norm_unit":score_inf, "final_fisher_step_inf_norm_unit":step_inf,
        "coarse_best_objective":float(best["objective"]),
        "coarse_second_best_gap":float(sols[1]["objective"]-best["objective"]),
        "evaluations_all_starts":int(sum(s["evaluations"] for s in sols)),
        "estimated_parameters_json":json.dumps([float(x) for x in refined["params"]]),
    }


def one_task(args):
    cfg, condition, rep, strategy, method_id = args
    angles = np.linspace(cfg["angle_min_deg"], cfg["angle_max_deg"], cfg["angle_bins"])
    gold = ROOT/cfg["material"]
    q = cfg["baseline_index"]+cfg["sample_delta_index"]
    d = cfg["gold_thickness_nm"]
    delta = condition["mismatch_nm"]
    tau = condition["target_angle_offset_deg"]
    aware = bool(condition["angle_aware"])
    allocation = cfg["allocations"]["angle_aware" if aware else "stable_angle"][strategy]
    stream = [cfg["seed"], condition["condition_id"], rep, method_id]
    rng = np.random.default_rng(np.random.SeedSequence(stream))
    if strategy == "T1":
        y = rng.poisson(truth_vector(cfg, angles, gold, gold, cfg["total_incident_photons"],
                                     [allocation[0]], q, d, delta, tau))
        fit = fit_one(y, "T1", cfg, angles, gold, [allocation[0]], fit_target_angle=aware)
    elif strategy == "T2":
        y = rng.poisson(truth_vector(cfg, angles, gold, gold, cfg["total_incident_photons"],
                                     allocation, q, d, delta, tau, t2=True))
        fit = fit_one(y, "T2", cfg, angles, gold, allocation, fit_target_angle=aware)
    elif strategy == "T_R":
        fref = allocation[1]
        y = rng.poisson(tr_truth(cfg, angles, gold, fref, d+delta, q, tau))
        fit = fit_tr(y, cfg, angles, gold, fref)
    else:
        raise ValueError(strategy)
    return {
        "task_id":f"c{condition['condition_id']}:r{rep}:{strategy}",
        "scenario_id":condition["condition_id"], "scenario":condition["scenario"],
        "replicate_id":rep, "strategy":strategy, "fit_stage":"final_refined",
        "config_hash":cfg["_config_hash"], "noise_stream_id":":".join(str(x) for x in stream),
        "true_mismatch_nm":delta, "true_target_angle_offset_deg":tau,
        "angle_aware":int(aware or strategy == "T_R"),
        "source_reference_fraction":allocation[0],
        "target_reference_fraction":allocation[1], "target_sample_fraction":allocation[2],
        "truth_q":q, **fit, "q_error_riu":fit["estimated_q"]-q,
    }


def summarize(rows):
    groups = {}
    for r in rows:
        key = (r["scenario_id"], r["strategy"])
        groups.setdefault(key, []).append(r)
    out = []
    for key in sorted(groups):
        group = groups[key]
        err = np.asarray([r["q_error_riu"] for r in group]); sq = err**2; n = len(err)
        out.append({
            "scenario_id":key[0], "scenario":group[0]["scenario"], "strategy":key[1], "n":n,
            "true_mismatch_nm":group[0]["true_mismatch_nm"],
            "true_target_angle_offset_deg":group[0]["true_target_angle_offset_deg"],
            "bias_riu":float(np.mean(err)), "bias_mcse_riu":float(np.std(err,ddof=1)/math.sqrt(n)),
            "variance_riu2":float(np.var(err,ddof=1)), "mse_riu2":float(np.mean(sq)),
            "mse_mcse_riu2":float(np.std(sq,ddof=1)/math.sqrt(n)),
            "rmse_riu":float(np.sqrt(np.mean(sq))), "mae_riu":float(np.mean(np.abs(err))),
            "p95_absolute_error_riu":float(np.quantile(np.abs(err),0.95)),
            "optimizer_failure_rate":float(np.mean([r["optimizer_failure"] for r in group])),
            "target_boundary_rate":float(np.mean([r["target_boundary"] for r in group])),
            "multi_solution_rate":float(np.mean([r["multi_solution"] for r in group])),
            "median_final_fisher_step_inf_norm_unit":float(np.median([r["final_fisher_step_inf_norm_unit"] for r in group])),
            "max_final_fisher_step_inf_norm_unit":float(np.max([r["final_fisher_step_inf_norm_unit"] for r in group])),
        })
    return out


def independent_differences(summaries):
    by = {(r["scenario_id"],r["strategy"]):r for r in summaries}
    out = []
    for cid in sorted({r["scenario_id"] for r in summaries}):
        for a, b in (("T1","T_R"),("T2","T_R"),("T1","T2")):
            ra, rb = by[(cid,a)], by[(cid,b)]
            mean = ra["mse_riu2"]-rb["mse_riu2"]
            se = math.sqrt(ra["mse_mcse_riu2"]**2+rb["mse_mcse_riu2"]**2)
            out.append({
                "scenario_id":cid, "scenario":ra["scenario"], "strategy_a":a, "strategy_b":b,
                "comparison_design":"independent protocol-specific Poisson streams; not paired",
                "mean_mse_difference_a_minus_b_riu2":mean, "mcse_riu2":se,
                "mc_interval_low_riu2":mean-1.96*se, "mc_interval_high_riu2":mean+1.96*se,
            })
    return out


def write_csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def main(config_path):
    cfg_path = Path(config_path).resolve()
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    cfg["_config_hash"] = digest(cfg_path)
    out = ROOT/"outputs"/"target_recalibration_confirmation"/cfg["run_id"]
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True)
    partial = out/"replicate_results.partial.csv"
    tasks = []
    for condition in cfg["conditions"]:
        for rep in range(cfg["replicates_per_condition"]):
            for method_id, strategy in enumerate(("T1","T2","T_R")):
                tasks.append((cfg, condition, rep, strategy, method_id))
    started = time.time(); rows = []; writer = None
    with open(partial,"w",newline="",encoding="utf-8") as f, ProcessPoolExecutor(max_workers=cfg["worker_processes"]) as pool:
        futures = [pool.submit(one_task, x) for x in tasks]
        for done, future in enumerate(as_completed(futures), start=1):
            row = future.result()
            if writer is None:
                writer = csv.DictWriter(f, fieldnames=list(row)); writer.writeheader()
            writer.writerow(row); f.flush(); rows.append(row)
            if done % 100 == 0 or done == len(futures):
                print(f"completed_tasks={done}/{len(futures)} records={len(rows)}", flush=True)
    rows.sort(key=lambda r:(r["scenario_id"],r["strategy"],r["replicate_id"]))
    write_csv(out/"replicate_results.csv",rows)
    summaries = summarize(rows); write_csv(out/"condition_summary.csv",summaries)
    differences = independent_differences(summaries); write_csv(out/"independent_risk_differences.csv",differences)
    summary = {
        "status":"FROZEN_TARGET_RECALIBRATION_CONFIRMATION_COMPLETE", "run_id":cfg["run_id"],
        "seed":cfg["seed"], "replicates_per_condition":cfg["replicates_per_condition"],
        "worker_processes":cfg["worker_processes"], "elapsed_seconds":time.time()-started,
        "config_sha256":cfg["_config_hash"], "script_sha256":digest(Path(__file__)),
        "material_sha256":digest(ROOT/cfg["material"]), "task_count":len(tasks),
        "record_count":len(rows), "allocation_selection_source":cfg["allocation_selection_source"],
        "allocation_selection_rule":cfg["allocation_selection_rule"],
        "comparison_note":"T1, T2, and T_R have different exposure protocols and independent Poisson streams. Their MSE differences use independent Monte Carlo standard errors and are not presented as paired differences.",
        "conditions":summaries, "differences":differences,
    }
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_target_recalibration_confirmation.py CONFIG")
    main(sys.argv[1])
