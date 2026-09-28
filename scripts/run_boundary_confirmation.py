#!/usr/bin/env python3
"""Frozen internal Monte Carlo confirmation for CalSPR V2.1."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"scripts"))
from calspr.inversion import bounded_pattern_fit, fisher_scoring_refine, poisson_deviance_half
from run_boundary_scan import bounds_for, model_vector, truth_vector


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def final_diagnostics(y, fn, params, bounds, h=1e-5):
    b = np.asarray(bounds, dtype=float)
    lower, width = b[:, 0], b[:, 1]-b[:, 0]
    u = np.clip((np.asarray(params)-lower)/width, 0.0, 1.0)
    mu = fn(params)
    jac = np.empty((len(mu), len(u)))
    for j in range(len(u)):
        plus, minus = u.copy(), u.copy()
        plus[j] = min(1.0, plus[j]+h)
        minus[j] = max(0.0, minus[j]-h)
        jac[:, j] = (fn(lower+plus*width)-fn(lower+minus*width))/(plus[j]-minus[j])
    score = jac.T@((np.asarray(y)-mu)/mu)
    fisher = jac.T@(jac/mu[:, None])
    try:
        step = np.linalg.solve(fisher+1e-12*np.diag(np.maximum(np.diag(fisher), 1.0)), score)
        step_inf = float(np.max(np.abs(step)))
    except np.linalg.LinAlgError:
        step_inf = float("inf")
    return float(np.max(np.abs(score))), step_inf


def fit_one(y, strategy, cfg, angles, gold, allocation, known_delta=0.0, fit_target_angle=False):
    bounds = bounds_for(strategy, cfg, fit_target_angle)
    fn = lambda p: model_vector(strategy, p, cfg, angles, gold,
                                cfg["total_incident_photons"], allocation,
                                known_delta=known_delta, fit_target_angle=fit_target_angle)
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
        "estimated_q": target,
        "objective_half_deviance": float(refined["objective"]),
        "optimizer_failure": int(not np.isfinite(refined["objective"])),
        "target_boundary": int(target-lo <= bf*(hi-lo) or hi-target <= bf*(hi-lo)),
        "multi_solution": int(multi),
        "refine_iterations": refined["iterations"],
        "refine_accepted_steps": refined["accepted_steps"],
        "final_score_inf_norm_unit": score_inf,
        "final_fisher_step_inf_norm_unit": step_inf,
        "coarse_best_objective": float(best["objective"]),
        "coarse_second_best_gap": float(sols[1]["objective"]-best["objective"]),
        "evaluations_all_starts": int(sum(s["evaluations"] for s in sols)),
        "estimated_parameters_json": json.dumps([float(x) for x in refined["params"]])
    }


def pair_task(args):
    cfg, condition_id, rep, delta, tau, angle_scenario = args
    angles = np.linspace(cfg["angle_min_deg"], cfg["angle_max_deg"], cfg["angle_bins"])
    gold = ROOT/cfg["material"]
    q = cfg["baseline_index"]+cfg["sample_delta_index"]
    allocation = [cfg["source_reference_fraction"]]
    stream = [cfg["seed"], condition_id, rep, 0]
    rng = np.random.default_rng(np.random.SeedSequence(stream))
    mu = truth_vector(cfg, angles, gold, gold, cfg["total_incident_photons"], allocation,
                      q, cfg["gold_thickness_nm"], delta, tau)
    y = rng.poisson(mu)
    strategies = (("T0", "T0", False), ("T1", "T1", False), ("O", "O", False))
    if angle_scenario:
        strategies = (("T1", "T1", False), ("T1_angle", "T1", True))
    rows = []
    for label, base, aware in strategies:
        fit = fit_one(y, base, cfg, angles, gold, allocation,
                      known_delta=delta if base == "O" else 0.0,
                      fit_target_angle=aware)
        rows.append({
            "scenario": "target_angle_offset" if angle_scenario else "thickness_mismatch",
            "condition_id": condition_id, "strategy": label, "replicate": rep,
            "noise_stream_id": ":".join(str(x) for x in stream),
            "true_mismatch_nm": delta, "true_target_angle_offset_deg": tau,
            "truth_q": q, **fit, "q_error_riu": fit["estimated_q"]-q
        })
    return rows


def t2_task(args):
    cfg, condition_id, rep, delta, tau, angle_scenario = args
    angles = np.linspace(cfg["angle_min_deg"], cfg["angle_max_deg"], cfg["angle_bins"])
    gold = ROOT/cfg["material"]
    q = cfg["baseline_index"]+cfg["sample_delta_index"]
    allocation = cfg["t2_allocation"]
    stream = [cfg["seed"], condition_id, rep, 1]
    rng = np.random.default_rng(np.random.SeedSequence(stream))
    mu = truth_vector(cfg, angles, gold, gold, cfg["total_incident_photons"], allocation,
                      q, cfg["gold_thickness_nm"], delta, tau, t2=True)
    y = rng.poisson(mu)
    choices = (("T2", False),) if not angle_scenario else (("T2", False), ("T2_angle", True))
    rows = []
    for label, aware in choices:
        fit = fit_one(y, "T2", cfg, angles, gold, allocation, fit_target_angle=aware)
        rows.append({
            "scenario": "target_angle_offset" if angle_scenario else "thickness_mismatch",
            "condition_id": condition_id, "strategy": label, "replicate": rep,
            "noise_stream_id": ":".join(str(x) for x in stream),
            "true_mismatch_nm": delta, "true_target_angle_offset_deg": tau,
            "truth_q": q, **fit, "q_error_riu": fit["estimated_q"]-q
        })
    return rows


def summarize(rows):
    groups = {}
    for r in rows:
        key = (r["scenario"], r["true_mismatch_nm"], r["true_target_angle_offset_deg"], r["strategy"])
        groups.setdefault(key, []).append(r)
    out = []
    for key, group in sorted(groups.items(), key=lambda z: tuple(str(x) for x in z[0])):
        err = np.array([r["q_error_riu"] for r in group])
        sq = err**2
        n = len(err)
        mse = float(np.mean(sq))
        out.append({
            "scenario": key[0], "true_mismatch_nm": key[1],
            "true_target_angle_offset_deg": key[2], "strategy": key[3], "n": n,
            "bias_riu": float(np.mean(err)), "bias_mcse_riu": float(np.std(err, ddof=1)/math.sqrt(n)),
            "variance_riu2": float(np.var(err, ddof=1)), "mse_riu2": mse,
            "mse_mcse_riu2": float(np.std(sq, ddof=1)/math.sqrt(n)),
            "rmse_riu": math.sqrt(mse), "mae_riu": float(np.mean(np.abs(err))),
            "p95_absolute_error_riu": float(np.quantile(np.abs(err), 0.95)),
            "optimizer_failure_rate": float(np.mean([r["optimizer_failure"] for r in group])),
            "target_boundary_rate": float(np.mean([r["target_boundary"] for r in group])),
            "multi_solution_rate": float(np.mean([r["multi_solution"] for r in group])),
            "median_final_fisher_step_inf_norm_unit": float(np.median([r["final_fisher_step_inf_norm_unit"] for r in group])),
            "max_final_fisher_step_inf_norm_unit": float(np.max([r["final_fisher_step_inf_norm_unit"] for r in group]))
        })
    return out


def main(config_path):
    cfg_path = Path(config_path).resolve()
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    out = ROOT/"outputs"/"boundary_confirmation"/cfg["run_id"]
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True)
    partial = out/"replicate_results.partial.csv"
    conditions = []
    cid = 0
    for delta in cfg["mismatch_nm"]:
        conditions.append((cid, delta, 0.0, False)); cid += 1
    for tau in cfg["target_angle_offset_deg"]:
        conditions.append((cid, 0.0, tau, True)); cid += 1
    tasks = []
    for cond in conditions:
        for rep in range(cfg["replicates_per_condition"]):
            args = (cfg, cond[0], rep, cond[1], cond[2], cond[3])
            tasks.append((pair_task, args))
            tasks.append((t2_task, args))
    started = time.time()
    rows = []
    writer = None
    with open(partial, "w", newline="", encoding="utf-8") as f, ProcessPoolExecutor(max_workers=cfg["worker_processes"]) as pool:
        futures = [pool.submit(fn, args) for fn, args in tasks]
        for done, future in enumerate(as_completed(futures), start=1):
            batch = future.result()
            for row in batch:
                if writer is None:
                    writer = csv.DictWriter(f, fieldnames=list(row)); writer.writeheader()
                writer.writerow(row); rows.append(row)
            f.flush()
            if done % 100 == 0 or done == len(futures):
                print(f"completed_tasks={done}/{len(futures)} records={len(rows)}", flush=True)
    rows.sort(key=lambda r: (r["condition_id"], r["strategy"], r["replicate"]))
    final = out/"replicate_results.csv"
    with open(final, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    partial.unlink()
    summaries = summarize(rows)
    with open(out/"condition_summary.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(summaries[0])); w.writeheader(); w.writerows(summaries)
    summary = {
        "status": "FROZEN_INTERNAL_CONFIRMATION_COMPLETE", "run_id": cfg["run_id"],
        "seed": cfg["seed"], "replicates_per_condition": cfg["replicates_per_condition"],
        "worker_processes": cfg["worker_processes"], "elapsed_seconds": time.time()-started,
        "config_sha256": digest(cfg_path), "material_sha256": digest(ROOT/cfg["material"]),
        "selection_source_run": cfg["selection_source_run"], "selection_rule": cfg["selection_rule"],
        "predeclared_interpretation": cfg["predeclared_interpretation"],
        "conditions": summaries
    }
    (out/"summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_boundary_confirmation.py CONFIG")
    main(sys.argv[1])
