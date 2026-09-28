#!/usr/bin/env python3
"""CalSPR V2.2 target-only recalibration strong-control screen."""
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
from calspr.inversion import bounded_pattern_fit, fisher_scoring_refine, poisson_deviance_half
from run_boundary_scan import (effective_info, model_vector, truth_vector,
                               t2_local_row, whitened_columns, local_row)


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def exposure(cfg, angles, gold, photons, d_nm, q, angle_zero=0.0, log_gain=0.0):
    return means(angles, cfg["wavelength_nm"], d_nm, q, gold, photons,
                 cfg["background_fraction"], angle_zero, log_gain)


def tr_model(p, cfg, angles, gold, target_reference_fraction):
    q, target_d, target_angle, log_gain = p
    fref = float(target_reference_fraction)
    budget = cfg["total_incident_photons"]
    return np.concatenate([
        exposure(cfg, angles, gold, budget*fref, target_d, cfg["baseline_index"], target_angle, log_gain),
        exposure(cfg, angles, gold, budget*(1.0-fref), target_d, q, target_angle, log_gain),
    ])


def tr_truth(cfg, angles, gold, target_reference_fraction, target_d, q, target_angle):
    return tr_model(np.array([q, target_d, target_angle, 0.0]), cfg, angles, gold,
                    target_reference_fraction)


def tr_bounds(cfg):
    b = cfg["bounds"]
    return [b["q"], b["source_thickness_nm"], b["angle_zero_deg"], b["log_gain"]]


def tr_local_row(cfg, angles, gold, target_reference_fraction):
    q = cfg["baseline_index"] + cfg["sample_delta_index"]
    p = np.array([q, cfg["gold_thickness_nm"], 0.0, 0.0])
    hs = cfg["finite_difference_steps"]
    steps = [hs["refractive_index"], hs["thickness_nm"], hs["angle_zero_deg"], hs["log_gain"]]
    fn = lambda z: tr_model(z, cfg, angles, gold, target_reference_fraction)
    _, cols = whitened_columns(fn, p, steps)
    info, _, rank, retained = effective_info(cols)
    return {
        "strategy": "T_R", "source_reference_fraction": 0.0,
        "target_reference_fraction": target_reference_fraction,
        "target_sample_fraction": 1.0-target_reference_fraction,
        "angle_model": "target_common_angle_fitted",
        "effective_information": info, "local_sd_riu": math.sqrt(1.0/info),
        "target_derivative_norm_retained": retained,
        "target_information_fraction_retained": retained**2,
        "nuisance_rank": rank,
    }


def fit_generic(y, fn, bounds, settings):
    sols = bounded_pattern_fit(lambda p: poisson_deviance_half(y, fn(p)), bounds, settings)
    best = sols[0]
    refined = fisher_scoring_refine(
        y, fn, best["params"], bounds,
        max_iter=settings["refine_max_iterations"],
        derivative_step=settings["refine_unit_derivative_step"],
    )
    return refined, float(sols[1]["objective"]-best["objective"])


def write_csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def normalize_local_rows(cfg, angles, gold):
    rows = []
    budget = cfg["total_incident_photons"]
    d = cfg["gold_thickness_nm"]
    dq = cfg["sample_delta_index"]
    for fs in cfg["t1_source_reference_fractions"]:
        r = local_row(cfg, angles, gold, budget, d, dq, fs)
        for mode, info_key, sd_key, retained_key, frac_key, rank_key in (
            ("stable_angle", "T1_effective_information", "T1_local_sd_riu", "T1_target_derivative_norm_retained", "T1_target_information_fraction_retained", "T1_nuisance_rank"),
            ("angle_aware", "T1_angle_aware_effective_information", "T1_angle_aware_local_sd_riu", "T1_angle_aware_target_derivative_norm_retained", "T1_angle_aware_target_information_fraction_retained", "T1_angle_aware_nuisance_rank"),
        ):
            rows.append({"strategy":"T1", "source_reference_fraction":fs,
                         "target_reference_fraction":0.0, "target_sample_fraction":1.0-fs,
                         "angle_model":mode, "effective_information":r[info_key],
                         "local_sd_riu":r[sd_key], "target_derivative_norm_retained":r[retained_key],
                         "target_information_fraction_retained":r[frac_key], "nuisance_rank":r[rank_key]})
    for alloc in cfg["t2_allocations"]:
        r = t2_local_row(cfg, angles, gold, budget, d, dq, alloc)
        for mode, info_key, sd_key, retained_key, frac_key, rank_key in (
            ("stable_angle", "T2_effective_information", "T2_local_sd_riu", "T2_target_derivative_norm_retained", "T2_target_information_fraction_retained", "T2_nuisance_rank"),
            ("angle_aware", "T2_angle_aware_effective_information", "T2_angle_aware_local_sd_riu", "T2_angle_aware_target_derivative_norm_retained", "T2_angle_aware_target_information_fraction_retained", "T2_angle_aware_nuisance_rank"),
        ):
            rows.append({"strategy":"T2", "source_reference_fraction":alloc[0],
                         "target_reference_fraction":alloc[1], "target_sample_fraction":alloc[2],
                         "angle_model":mode, "effective_information":r[info_key],
                         "local_sd_riu":r[sd_key], "target_derivative_norm_retained":r[retained_key],
                         "target_information_fraction_retained":r[frac_key], "nuisance_rank":r[rank_key]})
    for fref in cfg["tr_target_reference_fractions"]:
        rows.append(tr_local_row(cfg, angles, gold, fref))
    return rows


def select_allocations(rows):
    selected = {}
    for strategy in ("T1", "T2"):
        for mode in ("stable_angle", "angle_aware"):
            candidates = [r for r in rows if r["strategy"] == strategy and r["angle_model"] == mode]
            selected[f"{strategy}_{mode}"] = min(candidates, key=lambda r: r["local_sd_riu"])
    tr = min((r for r in rows if r["strategy"] == "T_R"), key=lambda r: r["local_sd_riu"])
    selected["T_R_target_common_angle_fitted"] = tr
    return selected


def deterministic_rows(cfg, angles, gold, selected):
    q = cfg["baseline_index"] + cfg["sample_delta_index"]
    d = cfg["gold_thickness_nm"]
    budget = cfg["total_incident_photons"]
    settings = cfg["fit"]
    out = []
    for state in cfg["representative_states"]:
        delta = state["mismatch_nm"]
        tau = state["target_angle_offset_deg"]
        aware = bool(state["angle_aware"])
        mode = "angle_aware" if aware else "stable_angle"
        for strategy in ("T1", "T2", "T_R"):
            key = f"{strategy}_{mode}" if strategy != "T_R" else "T_R_target_common_angle_fitted"
            sel = selected[key]
            allocation = [sel["source_reference_fraction"], sel["target_reference_fraction"], sel["target_sample_fraction"]]
            if strategy == "T1":
                alloc = [allocation[0]]
                y = truth_vector(cfg, angles, gold, gold, budget, alloc, q, d, delta, tau)
                if aware:
                    bounds = [cfg["bounds"]["q"], cfg["bounds"]["source_thickness_nm"], cfg["bounds"]["mismatch_nm"], cfg["bounds"]["target_angle_offset_deg"], cfg["bounds"]["angle_zero_deg"], cfg["bounds"]["log_gain"]]
                else:
                    bounds = [cfg["bounds"]["q"], cfg["bounds"]["source_thickness_nm"], cfg["bounds"]["mismatch_nm"], cfg["bounds"]["angle_zero_deg"], cfg["bounds"]["log_gain"]]
                fn = lambda p, a=alloc, aw=aware: model_vector("T1", p, cfg, angles, gold, budget, a, fit_target_angle=aw)
            elif strategy == "T2":
                y = truth_vector(cfg, angles, gold, gold, budget, allocation, q, d, delta, tau, t2=True)
                if aware:
                    bounds = [cfg["bounds"]["q"], cfg["bounds"]["source_thickness_nm"], cfg["bounds"]["mismatch_nm"], cfg["bounds"]["target_angle_offset_deg"], cfg["bounds"]["angle_zero_deg"], cfg["bounds"]["log_gain"]]
                else:
                    bounds = [cfg["bounds"]["q"], cfg["bounds"]["source_thickness_nm"], cfg["bounds"]["mismatch_nm"], cfg["bounds"]["angle_zero_deg"], cfg["bounds"]["log_gain"]]
                fn = lambda p, a=allocation, aw=aware: model_vector("T2", p, cfg, angles, gold, budget, a, fit_target_angle=aw)
            else:
                fref = allocation[1]
                y = tr_truth(cfg, angles, gold, fref, d+delta, q, tau)
                bounds = tr_bounds(cfg)
                fn = lambda p, fr=fref: tr_model(p, cfg, angles, gold, fr)
            fit, gap = fit_generic(y, fn, bounds, settings)
            est = float(fit["params"][0])
            out.append({
                "state_id":state["state_id"], "strategy":strategy,
                "true_mismatch_nm":delta, "true_target_angle_offset_deg":tau,
                "angle_aware":int(aware or strategy == "T_R"),
                "source_reference_fraction":allocation[0],
                "target_reference_fraction":allocation[1],
                "target_sample_fraction":allocation[2],
                "truth_q":q, "estimated_q":est, "q_bias_riu":est-q,
                "objective_half_deviance":float(fit["objective"]),
                "coarse_second_best_gap":gap,
                "estimated_parameters_json":json.dumps([float(x) for x in fit["params"]]),
            })
    return out


def main(config_path):
    cfg_path = Path(config_path).resolve()
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    out = ROOT / "outputs" / "target_recalibration_scan" / cfg["run_id"]
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True)
    started = time.time()
    angles = np.linspace(cfg["angle_min_deg"], cfg["angle_max_deg"], cfg["angle_bins"])
    gold = ROOT / cfg["material"]
    rows = normalize_local_rows(cfg, angles, gold)
    write_csv(out / "allocation_scan.csv", rows)
    selected = select_allocations(rows)
    deterministic = deterministic_rows(cfg, angles, gold, selected)
    write_csv(out / "deterministic_representative.csv", deterministic)
    selection = {
        "selection_rule": cfg["selection_rule"],
        "selected": selected,
        "noisy_outcomes_used_for_selection": False,
    }
    (out / "selection.json").write_text(json.dumps(selection, indent=2)+"\n", encoding="utf-8")
    summary = {
        "status":"EXPLORATORY_STRONG_CONTROL_SCREEN_COMPLETE",
        "run_id":cfg["run_id"], "seed":cfg["seed"],
        "elapsed_seconds":time.time()-started,
        "config_sha256":digest(cfg_path), "script_sha256":digest(Path(__file__)),
        "material_sha256":digest(gold), "allocation_rows":len(rows),
        "deterministic_rows":len(deterministic), "selected":selected,
        "interpretation_limit":cfg["interpretation_limit"],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_target_recalibration_scan.py CONFIG")
    main(sys.argv[1])
