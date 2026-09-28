#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from calspr.information import means
from calspr.inversion import bounded_pattern_fit, poisson_deviance_half


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def model_vector(regime, p, cfg, angles, gold):
    w, total, bg, q0 = cfg["wavelength_nm"], cfg["total_incident_photons"], cfg["background_fraction"], cfg["baseline_index"]
    if regime == "absolute":
        q, d, a0, lg = p
        return means(angles, w, d, q, gold, total, bg, a0, lg)
    if regime == "same_device_differential":
        delta, d, a0, lg = p
        return np.concatenate([means(angles, w, d, q0, gold, total/2, bg, a0, lg),
                               means(angles, w, d, q0+delta, gold, total/2, bg, a0, lg)])
    if regime == "transfer_naive":
        q, d, a0, lg = p
        return np.concatenate([means(angles, w, d, q0, gold, total/2, bg, a0, lg),
                               means(angles, w, d, q, gold, total/2, bg, a0, lg)])
    if regime == "transfer_aware":
        q, d, mismatch, a0, lg = p
        return np.concatenate([means(angles, w, d, q0, gold, total/2, bg, a0, lg),
                               means(angles, w, d+mismatch, q, gold, total/2, bg, a0, lg)])
    raise ValueError(regime)


def truth_vector(regime, mismatch, cfg, angles, gold):
    q0, dq, d = cfg["baseline_index"], cfg["sample_delta_index"], cfg["gold_thickness_nm"]
    w, total, bg = cfg["wavelength_nm"], cfg["total_incident_photons"], cfg["background_fraction"]
    if regime == "absolute":
        return means(angles, w, d+mismatch, q0+dq, gold, total, bg)
    if regime == "same_device_differential":
        return np.concatenate([means(angles, w, d, q0, gold, total/2, bg),
                               means(angles, w, d, q0+dq, gold, total/2, bg)])
    return np.concatenate([means(angles, w, d, q0, gold, total/2, bg),
                           means(angles, w, d+mismatch, q0+dq, gold, total/2, bg)])


def target_truth(regime, cfg):
    return cfg["sample_delta_index"] if regime == "same_device_differential" else cfg["baseline_index"]+cfg["sample_delta_index"]


def summarize(group):
    errors = np.array([r["target_error_riu"] for r in group])
    return {
        "n": len(group), "bias_riu": float(np.mean(errors)), "rmse_riu": float(np.sqrt(np.mean(errors**2))),
        "mae_riu": float(np.mean(np.abs(errors))), "p95_absolute_error_riu": float(np.quantile(np.abs(errors), 0.95)),
        "optimizer_failure_rate": float(np.mean([r["optimizer_failure"] for r in group])),
        "target_boundary_rate": float(np.mean([r["target_boundary"] for r in group])),
        "multi_solution_rate": float(np.mean([r["multi_solution"] for r in group])),
        "median_evaluations_all_starts": float(np.median([r["evaluations_all_starts"] for r in group]))
    }


def main(config_path):
    cfg_path = Path(config_path).resolve()
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    out = ROOT / "outputs" / "p3_inversion" / cfg["run_id"]
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True)
    gold = ROOT / cfg["material_sources"]["gold_raw"]
    angles = np.linspace(cfg["angle_min_deg"], cfg["angle_max_deg"], cfg["angle_bins"])
    rng = np.random.default_rng(cfg["seed"])
    regimes = ["absolute", "same_device_differential", "transfer_naive", "transfer_aware"]
    conditions = [("absolute", 0.0), ("same_device_differential", 0.0)]
    conditions += [(regime, m) for regime in ("transfer_naive", "transfer_aware") for m in cfg["transfer_mismatch_nm"]]
    records = []
    start = time.time()
    for regime, mismatch in conditions:
        truth_mu = truth_vector(regime, mismatch, cfg, angles, gold)
        truth_target = target_truth(regime, cfg)
        bounds = cfg["bounds"][regime]
        for rep in range(cfg["replicates_per_condition"]):
            y = rng.poisson(truth_mu)
            solutions = bounded_pattern_fit(lambda p: poisson_deviance_half(y, model_vector(regime, p, cfg, angles, gold)), bounds, cfg["fit"])
            best = solutions[0]
            target = float(best["params"][0])
            delta_obj = cfg["fit"]["multi_solution_deviance_delta"]
            sep = cfg["fit"]["multi_solution_target_separation_riu"]
            multi = any((s["objective"]-best["objective"] <= delta_obj and abs(float(s["params"][0])-target) >= sep) for s in solutions[1:])
            lo, hi = bounds[0]
            bfrac = cfg["fit"]["target_boundary_fraction"]
            boundary = (target-lo <= bfrac*(hi-lo)) or (hi-target <= bfrac*(hi-lo))
            failure = (not np.isfinite(best["objective"])) or best["final_step"] > cfg["fit"]["normalized_step_tolerance"]
            records.append({
                "regime": regime, "true_transfer_mismatch_nm": mismatch, "replicate": rep,
                "truth_target": truth_target, "estimated_target": target, "target_error_riu": target-truth_target,
                "objective_half_deviance": best["objective"], "optimizer_failure": int(failure),
                "target_boundary": int(boundary), "multi_solution": int(multi),
                "evaluations_all_starts": sum(s["evaluations"] for s in solutions),
                "best_parameters_json": json.dumps([float(x) for x in best["params"]])
            })
        print(f"completed regime={regime} mismatch_nm={mismatch}", flush=True)

    with open(out/"replicate_results.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0])); writer.writeheader(); writer.writerows(records)
    grouped = []
    for regime, mismatch in conditions:
        g = [r for r in records if r["regime"] == regime and r["true_transfer_mismatch_nm"] == mismatch]
        grouped.append({"regime": regime, "true_transfer_mismatch_nm": mismatch, **summarize(g)})
    with open(out/"condition_summary.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(grouped[0])); writer.writeheader(); writer.writerows(grouped)
    summary = {
        "status": "EXPLORATORY_COMPLETE", "run_id": cfg["run_id"], "seed": cfg["seed"],
        "elapsed_seconds": time.time()-start, "config_sha256": sha256(cfg_path), "gold_source_sha256": sha256(gold),
        "virtual_cluster_unit": "one Poisson replicate for one device or device-pair condition",
        "initialization": "deterministic Halton library over full bounds; independent of truth",
        "conditions": grouped,
        "interpretation_limit": "scenario-based computational inversion; material and drift distributions are not experimentally validated"
    }
    (out/"summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_p3_inversion.py CONFIG")
    main(sys.argv[1])

