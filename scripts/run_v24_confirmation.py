#!/usr/bin/env python3
"""Paired T0/T1 confirmation of nine precomputed V2.4 threshold anchors."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from calspr.inversion import bounded_pattern_fit, fisher_scoring_refine, poisson_deviance_half
from calspr.v24_design import ResponseScenario, angle_centres, gold_reflectance, poisson_mean


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def exposure(cfg, angles, gold, scenario, photons, d, q, a0=0.0, lg=0.0):
    r = gold_reflectance(angles+a0, cfg["wavelength_nm"], d, q, gold, scenario)
    return poisson_mean(r, photons, cfg["background_fraction"], lg)


def truth(cfg, angles, gold, scenario, anchor, delta):
    q0 = cfg["baseline_index"]
    q = q0 + cfg["sample_delta_index"]
    d = anchor["thickness_nm"]
    budget = anchor["total_incident_photons"]
    return np.concatenate([
        exposure(cfg, angles, gold, scenario, budget*0.5, d, q0),
        exposure(cfg, angles, gold, scenario, budget*0.5, d+delta, q),
    ])


def model(cfg, angles, gold, scenario, anchor, strategy, p):
    q0 = cfg["baseline_index"]
    budget = anchor["total_incident_photons"]
    if strategy == "T0":
        q, d, a0, lg = p
        delta = 0.0
    else:
        q, d, delta, a0, lg = p
    return np.concatenate([
        exposure(cfg, angles, gold, scenario, budget*0.5, d, q0, a0, lg),
        exposure(cfg, angles, gold, scenario, budget*0.5, d+delta, q, a0, lg),
    ])


def bounds_for(cfg, strategy):
    b = cfg["bounds"]
    if strategy == "T0":
        return [b["q"], b["source_thickness_nm"], b["angle_zero_deg"], b["log_gain"]]
    return [b["q"], b["source_thickness_nm"], b["mismatch_nm"],
            b["angle_zero_deg"], b["log_gain"]]


def final_diagnostics(y, fn, params, bounds, h=1e-5):
    bounds = np.asarray(bounds, dtype=float)
    lower, width = bounds[:, 0], bounds[:, 1]-bounds[:, 0]
    u = np.clip((np.asarray(params)-lower)/width, 0, 1)
    mu = fn(params)
    jac = np.empty((len(mu), len(u)))
    for j in range(len(u)):
        plus, minus = u.copy(), u.copy()
        plus[j] = min(1.0, plus[j]+h)
        minus[j] = max(0.0, minus[j]-h)
        jac[:, j] = (fn(lower+plus*width)-fn(lower+minus*width))/(plus[j]-minus[j])
    score = jac.T @ ((np.asarray(y)-mu)/mu)
    fisher = jac.T @ (jac/mu[:, None])
    try:
        step = np.linalg.solve(fisher+1e-12*np.diag(np.maximum(np.diag(fisher), 1.0)), score)
        step_inf = float(np.max(np.abs(step)))
    except np.linalg.LinAlgError:
        step_inf = math.inf
    score_scaled = float(np.max(np.abs(score)/np.sqrt(np.maximum(np.diag(fisher), 1e-300))))
    return score_scaled, step_inf


def fit_one(y, cfg, angles, gold, scenario, anchor, strategy):
    bounds = bounds_for(cfg, strategy)
    fn = lambda p: model(cfg, angles, gold, scenario, anchor, strategy, p)
    sols = bounded_pattern_fit(lambda p: poisson_deviance_half(y, fn(p)), bounds, cfg["fit"])
    best = sols[0]
    refined = fisher_scoring_refine(
        y, fn, best["params"], bounds,
        max_iter=cfg["fit"]["refine_max_iterations"],
        derivative_step=cfg["fit"]["refine_unit_derivative_step"],
    )
    score, step = final_diagnostics(y, fn, refined["params"], bounds,
                                    cfg["fit"]["refine_unit_derivative_step"])
    q_fraction = (refined["params"][0]-bounds[0][0])/(bounds[0][1]-bounds[0][0])
    multi = any(
        s["objective"]-best["objective"] <= cfg["fit"]["multi_solution_deviance_delta"] and
        abs(float(s["params"][0])-float(best["params"][0])) >= cfg["fit"]["multi_solution_target_separation_riu"]
        for s in sols[1:]
    )
    return {
        "estimated_q": float(refined["params"][0]),
        "estimated_parameters_json": json.dumps([float(x) for x in refined["params"]]),
        "objective_half_deviance": float(refined["objective"]),
        "optimizer_failure": int(not np.isfinite(refined["objective"]) or not np.isfinite(step)),
        "target_boundary": int(q_fraction <= cfg["fit"]["target_boundary_fraction"] or q_fraction >= 1-cfg["fit"]["target_boundary_fraction"]),
        "multi_solution": int(multi),
        "final_score_scaled_inf": score,
        "final_fisher_step_inf_norm_unit": step,
        "fit_iterations": refined["iterations"],
        "fit_accepted_steps": refined["accepted_steps"],
    }


def paired_task(args):
    cfg, anchor, rho, rep = args
    angles = angle_centres(cfg)
    gold = ROOT / cfg["material"]
    a = cfg["angular_response"]
    scenario = ResponseScenario(a["name"], a["kind"], a["width_deg"], a["quadrature_order"])
    delta = rho * anchor["delta_pred_nm"]
    anchor_index = next(i for i, x in enumerate(cfg["anchors"]) if x["anchor_id"] == anchor["anchor_id"])
    rho_index = cfg["relative_magnitudes"].index(rho)
    stream = [cfg["seed"], anchor_index, rho_index, rep]
    rng = np.random.default_rng(np.random.SeedSequence(stream))
    y = rng.poisson(truth(cfg, angles, gold, scenario, anchor, delta))
    qtrue = cfg["baseline_index"] + cfg["sample_delta_index"]
    out = []
    for strategy in ("T0", "T1"):
        fit = fit_one(y, cfg, angles, gold, scenario, anchor, strategy)
        err = fit["estimated_q"]-qtrue
        out.append({
            "task_id": f"{anchor['anchor_id']}|rho={rho:+g}|rep={rep:03d}",
            "anchor_id": anchor["anchor_id"], "thickness_nm": anchor["thickness_nm"],
            "total_incident_photons": anchor["total_incident_photons"],
            "rho_delta_over_pred": rho, "delta_pred_nm": anchor["delta_pred_nm"],
            "true_delta_nm": delta, "replicate": rep, "strategy": strategy,
            "noise_stream_json": json.dumps(stream), "true_q": qtrue,
            "error_riu": err, "squared_error_riu2": err*err,
            **fit,
        })
    return out


def summarize(rows, cfg):
    groups = defaultdict(list)
    for row in rows:
        groups[(row["anchor_id"], float(row["rho_delta_over_pred"]), row["strategy"])].append(row)
    summaries = []
    for (anchor_id, rho, strategy), group in sorted(groups.items()):
        err = np.asarray([r["error_riu"] for r in group])
        sq = err*err
        n = len(group)
        summaries.append({
            "summary_type": "strategy", "anchor_id": anchor_id,
            "rho_delta_over_pred": rho, "strategy": strategy, "n": n,
            "bias_riu": float(np.mean(err)),
            "bias_mcse_riu": float(np.std(err, ddof=1)/math.sqrt(n)),
            "mse_riu2": float(np.mean(sq)),
            "rmse_riu": float(np.sqrt(np.mean(sq))),
            "failure_rate": float(np.mean([r["optimizer_failure"] for r in group])),
            "target_boundary_rate": float(np.mean([r["target_boundary"] for r in group])),
            "multi_solution_rate": float(np.mean([r["multi_solution"] for r in group])),
        })
    by_task = defaultdict(dict)
    for row in rows:
        by_task[row["task_id"]][row["strategy"]] = row
    paired = defaultdict(list)
    for task_id, item in by_task.items():
        if set(item) != {"T0", "T1"}:
            raise RuntimeError(f"incomplete pair {task_id}")
        key = (item["T0"]["anchor_id"], float(item["T0"]["rho_delta_over_pred"]))
        paired[key].append(item["T0"]["squared_error_riu2"]-item["T1"]["squared_error_riu2"])
    anchor_map = {a["anchor_id"]: a for a in cfg["anchors"]}
    extension_triggered = False
    for (anchor_id, rho), values in sorted(paired.items()):
        values = np.asarray(values)
        n = len(values)
        mean = float(np.mean(values))
        mcse = float(np.std(values, ddof=1)/math.sqrt(n))
        low, high = mean-1.96*mcse, mean+1.96*mcse
        a = anchor_map[anchor_id]
        predicted = (a["v1_riu2"]-a["v0_riu2"])*(rho*rho-1.0)
        trigger = bool(abs(rho) in (0.5, 2.0) and low <= 0 <= high and
                       mcse > 0.25*abs(predicted))
        extension_triggered = extension_triggered or trigger
        summaries.append({
            "summary_type": "paired_D_T0_minus_T1", "anchor_id": anchor_id,
            "rho_delta_over_pred": rho, "strategy": "T0-minus-T1", "n": n,
            "mean_D_riu2": mean, "D_mcse_riu2": mcse,
            "D_mc_interval_low_riu2": low, "D_mc_interval_high_riu2": high,
            "local_predicted_D_riu2": predicted,
            "extension_rule_trigger": int(trigger),
        })
    return summaries, extension_triggered


def write_csv(path, rows):
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader(); writer.writerows(rows)


def main(config_path):
    started = time.time()
    cfg_path = Path(config_path).resolve()
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    source = ROOT / cfg["selection_source"]
    if digest(source) != cfg["selection_source_sha256"]:
        raise RuntimeError("threshold prediction source hash mismatch")
    out = ROOT / "outputs" / "v24" / cfg["run_id"]
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True)
    tasks = [
        (cfg, anchor, float(rho), rep)
        for anchor in cfg["anchors"]
        for rho in cfg["relative_magnitudes"]
        for rep in range(cfg["replicates_per_condition"])
    ]
    partial = out / "V24_REPLICATES.partial.csv"
    rows = []
    with open(partial, "w", newline="", encoding="utf-8") as f, ProcessPoolExecutor(max_workers=cfg["worker_processes"]) as pool:
        writer = None
        futures = [pool.submit(paired_task, task) for task in tasks]
        for done, future in enumerate(as_completed(futures), start=1):
            pair = future.result()
            if writer is None:
                writer = csv.DictWriter(f, fieldnames=list(pair[0]))
                writer.writeheader()
            writer.writerows(pair); f.flush(); rows.extend(pair)
            if done % 100 == 0:
                print(f"completed paired tasks {done}/{len(tasks)}", flush=True)
    rows.sort(key=lambda r: (r["anchor_id"], float(r["rho_delta_over_pred"]), r["replicate"], r["strategy"]))
    final_replicates = out / "V24_REPLICATES.csv"
    write_csv(final_replicates, rows)
    partial.unlink()
    summaries, extend = summarize(rows, cfg)
    write_csv(out / "V24_CONFIRMATION_SUMMARY.csv", summaries)
    expected = len(tasks)*2
    summary = {
        "status": "V24_CONFIRMATION_COMPLETE",
        "run_id": cfg["run_id"], "seed": cfg["seed"],
        "replicates_per_condition": cfg["replicates_per_condition"],
        "worker_processes": cfg["worker_processes"],
        "task_count": len(tasks), "record_count": len(rows), "expected_record_count": expected,
        "record_count_pass": len(rows) == expected,
        "optimizer_failure_count": int(sum(r["optimizer_failure"] for r in rows)),
        "target_boundary_count": int(sum(r["target_boundary"] for r in rows)),
        "multi_solution_count": int(sum(r["multi_solution"] for r in rows)),
        "extension_rule_triggered": bool(extend),
        "extension_rule": cfg["extension_rule"],
        "elapsed_seconds": time.time()-started,
        "peak_process_count": cfg["worker_processes"],
        "config_sha256": digest(cfg_path),
        "selection_source_sha256": digest(source),
        "replicates_sha256": digest(final_replicates),
        "interpretation_limit": cfg["interpretation_limit"],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_v24_confirmation.py CONFIG")
    main(sys.argv[1])
