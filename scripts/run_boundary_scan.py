#!/usr/bin/env python3
"""Exploratory calibration-transfer boundary analysis for CalSPR V2.1.

This run separates three products: local derivative geometry, deterministic
pseudo-true fits under finite mismatch, and representative material-model
mismatch.  It does not generate random Monte Carlo confirmation data.
"""
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


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def exposure(cfg, angles, gold, photons, d_nm, q, angle_zero=0.0, log_gain=0.0):
    return means(angles, cfg["wavelength_nm"], d_nm, q, gold, photons,
                 cfg["background_fraction"], angle_zero, log_gain)


def model_vector(strategy, p, cfg, angles, gold, budget, allocation, known_delta=0.0,
                 fit_target_angle=False):
    q0 = cfg["baseline_index"]
    if strategy in ("T0", "O"):
        q, d, a0, lg = p
        delta = known_delta if strategy == "O" else 0.0
        fs = allocation[0]
        return np.concatenate([
            exposure(cfg, angles, gold, budget*fs, d, q0, a0, lg),
            exposure(cfg, angles, gold, budget*(1-fs), d+delta, q, a0, lg)])
    if strategy == "T1":
        if fit_target_angle:
            q, d, delta, tau, a0, lg = p
        else:
            q, d, delta, a0, lg = p
            tau = 0.0
        fs = allocation[0]
        return np.concatenate([
            exposure(cfg, angles, gold, budget*fs, d, q0, a0, lg),
            exposure(cfg, angles, gold, budget*(1-fs), d+delta, q, a0+tau, lg)])
    if strategy == "T2":
        if fit_target_angle:
            q, d, delta, tau, a0, lg = p
        else:
            q, d, delta, a0, lg = p
            tau = 0.0
        fs, ftr, fts = allocation
        return np.concatenate([
            exposure(cfg, angles, gold, budget*fs, d, q0, a0, lg),
            exposure(cfg, angles, gold, budget*ftr, d+delta, q0, a0+tau, lg),
            exposure(cfg, angles, gold, budget*fts, d+delta, q, a0+tau, lg)])
    raise ValueError(strategy)


def truth_vector(cfg, angles, source_gold, target_gold, budget, allocation, q, d,
                 delta, tau=0.0, t2=False):
    q0 = cfg["baseline_index"]
    if t2:
        fs, ftr, fts = allocation
        return np.concatenate([
            exposure(cfg, angles, source_gold, budget*fs, d, q0),
            exposure(cfg, angles, target_gold, budget*ftr, d+delta, q0, tau),
            exposure(cfg, angles, target_gold, budget*fts, d+delta, q, tau)])
    fs = allocation[0]
    return np.concatenate([
        exposure(cfg, angles, source_gold, budget*fs, d, q0),
        exposure(cfg, angles, target_gold, budget*(1-fs), d+delta, q, tau)])


def derivative(fn, x, step):
    return (fn(x+step)-fn(x-step))/(2.0*step)


def residualize(target, nuisance, rtol=1e-10):
    target = np.asarray(target, dtype=float)
    nuisance = np.asarray(nuisance, dtype=float)
    if nuisance.size == 0:
        return target.copy(), 0, 1.0
    norms = np.linalg.norm(nuisance, axis=0)
    keep = norms > 0
    scaled = nuisance[:, keep] / norms[keep]
    u, s, _ = np.linalg.svd(scaled, full_matrices=False)
    tol = rtol * (s[0] if len(s) else 1.0)
    rank = int(np.sum(s > tol))
    resid = target-u[:, :rank]@(u[:, :rank].T@target)
    return resid, rank, float(np.linalg.norm(resid)/np.linalg.norm(target))


def whitened_columns(fn, p, steps):
    mu = fn(p)
    if np.any(mu <= 0):
        raise ValueError("non-positive Poisson mean")
    cols = []
    for j, h in enumerate(steps):
        cols.append(derivative(lambda x: fn(np.r_[p[:j], x, p[j+1:]]), p[j], h)/np.sqrt(mu))
    return mu, np.column_stack(cols)


def effective_info(cols):
    rq, rank, retained = residualize(cols[:, 0], cols[:, 1:])
    return float(rq@rq), rq, rank, retained


def local_row(cfg, angles, gold, budget, d, dq, ref_fraction):
    q = cfg["baseline_index"]+dq
    hs = cfg["finite_difference_steps"]
    p0 = np.array([q, d, 0.0, 0.0])
    steps0 = [hs["refractive_index"], hs["thickness_nm"], hs["angle_zero_deg"], hs["log_gain"]]
    f0 = lambda p: model_vector("T0", p, cfg, angles, gold, budget, [ref_fraction])
    mu0, c0 = whitened_columns(f0, p0, steps0)
    i0, rq, rank0, retain0 = effective_info(c0)
    wroot = 1.0/np.sqrt(mu0)
    rdelta_raw = derivative(
        lambda z: truth_vector(cfg, angles, gold, gold, budget, [ref_fraction], q, d, z),
        0.0, hs["thickness_nm"])*wroot
    rtau_raw = derivative(
        lambda z: truth_vector(cfg, angles, gold, gold, budget, [ref_fraction], q, d, 0.0, z),
        0.0, hs["angle_zero_deg"])*wroot
    nuisance0 = c0[:, 1:]
    rdelta, _, _ = residualize(rdelta_raw, nuisance0)
    rtau, _, _ = residualize(rtau_raw, nuisance0)
    kdelta = float(rq@rdelta/i0)
    ktau = float(rq@rtau/i0)

    p1 = np.array([q, d, 0.0, 0.0, 0.0])
    steps1 = [hs["refractive_index"], hs["thickness_nm"], hs["thickness_nm"],
              hs["angle_zero_deg"], hs["log_gain"]]
    f1 = lambda p: model_vector("T1", p, cfg, angles, gold, budget, [ref_fraction])
    _, c1 = whitened_columns(f1, p1, steps1)
    i1, rq1, rank1, retain1 = effective_info(c1)

    p1a = np.array([q, d, 0.0, 0.0, 0.0, 0.0])
    steps1a = [hs["refractive_index"], hs["thickness_nm"], hs["thickness_nm"],
               hs["angle_zero_deg"], hs["angle_zero_deg"], hs["log_gain"]]
    f1a = lambda p: model_vector("T1", p, cfg, angles, gold, budget, [ref_fraction], fit_target_angle=True)
    _, c1a = whitened_columns(f1a, p1a, steps1a)
    i1a, _, rank1a, retain1a = effective_info(c1a)

    c = float((rq@rdelta)/(np.linalg.norm(rq)*np.linalg.norm(rdelta))) if np.linalg.norm(rdelta) else 0.0
    predicted_i1 = i0*(1.0-c*c)
    var0 = 1.0/i0
    var1 = 1.0/i1
    dstar = math.sqrt(var1-var0)/abs(kdelta) if var1 >= var0 and kdelta != 0 else None
    return {
        "budget_photons": budget, "source_reference_fraction": ref_fraction,
        "sample_delta_index": dq, "nominal_thickness_nm": d,
        "T0_effective_information": i0, "T0_local_sd_riu": math.sqrt(var0),
        "T0_target_derivative_norm_retained": retain0,
        "T0_target_information_fraction_retained": retain0**2,
        "T0_nuisance_rank": rank0, "K_delta_riu_per_nm": kdelta,
        "K_target_angle_riu_per_deg": ktau,
        "T1_effective_information": i1, "T1_local_sd_riu": math.sqrt(var1),
        "T1_target_derivative_norm_retained": retain1,
        "T1_target_information_fraction_retained": retain1**2,
        "T1_nuisance_rank": rank1, "residual_angle_cosine_delta": c,
        "T1_information_from_angle_formula": predicted_i1,
        "angle_formula_relative_error": abs(predicted_i1-i1)/i1 if i1 else None,
        "T0_vs_T1_boundary_abs_delta_nm": dstar,
        "T1_angle_aware_effective_information": i1a,
        "T1_angle_aware_local_sd_riu": math.sqrt(1.0/i1a),
        "T1_angle_aware_target_derivative_norm_retained": retain1a,
        "T1_angle_aware_target_information_fraction_retained": retain1a**2,
        "T1_angle_aware_nuisance_rank": rank1a
    }


def t2_local_row(cfg, angles, gold, budget, d, dq, allocation):
    q = cfg["baseline_index"]+dq
    hs = cfg["finite_difference_steps"]
    p = np.array([q, d, 0.0, 0.0, 0.0])
    steps = [hs["refractive_index"], hs["thickness_nm"], hs["thickness_nm"],
             hs["angle_zero_deg"], hs["log_gain"]]
    fn = lambda z: model_vector("T2", z, cfg, angles, gold, budget, allocation)
    _, cols = whitened_columns(fn, p, steps)
    info, _, rank, retained = effective_info(cols)
    pa = np.array([q, d, 0.0, 0.0, 0.0, 0.0])
    stepsa = [hs["refractive_index"], hs["thickness_nm"], hs["thickness_nm"],
              hs["angle_zero_deg"], hs["angle_zero_deg"], hs["log_gain"]]
    fna = lambda z: model_vector("T2", z, cfg, angles, gold, budget, allocation, fit_target_angle=True)
    _, colsa = whitened_columns(fna, pa, stepsa)
    infoa, _, ranka, retaineda = effective_info(colsa)
    return {
        "budget_photons": budget, "source_reference_fraction": allocation[0],
        "target_reference_fraction": allocation[1], "target_sample_fraction": allocation[2],
        "sample_delta_index": dq, "nominal_thickness_nm": d,
        "T2_effective_information": info, "T2_local_sd_riu": math.sqrt(1.0/info),
        "T2_target_derivative_norm_retained": retained,
        "T2_target_information_fraction_retained": retained**2,
        "T2_nuisance_rank": rank,
        "T2_angle_aware_effective_information": infoa,
        "T2_angle_aware_local_sd_riu": math.sqrt(1.0/infoa),
        "T2_angle_aware_target_derivative_norm_retained": retaineda,
        "T2_angle_aware_target_information_fraction_retained": retaineda**2,
        "T2_angle_aware_nuisance_rank": ranka
    }


def bounds_for(strategy, cfg, fit_target_angle=False):
    b = cfg["bounds"]
    if strategy in ("T0", "O"):
        return [b["q"], b["source_thickness_nm"], b["angle_zero_deg"], b["log_gain"]]
    out = [b["q"], b["source_thickness_nm"], b["mismatch_nm"]]
    if fit_target_angle:
        out.append(b["target_angle_offset_deg"])
    return out+[b["angle_zero_deg"], b["log_gain"]]


def deterministic_fit(strategy, y, cfg, angles, fit_gold, budget, allocation,
                      known_delta=0.0, fit_target_angle=False):
    bounds = bounds_for(strategy, cfg, fit_target_angle)
    settings = cfg["fit"]
    fn = lambda p: model_vector(strategy, p, cfg, angles, fit_gold, budget, allocation,
                                known_delta=known_delta, fit_target_angle=fit_target_angle)
    sols = bounded_pattern_fit(lambda p: poisson_deviance_half(y, fn(p)), bounds, settings)
    refined = fisher_scoring_refine(y, fn, sols[0]["params"], bounds, max_iter=60, derivative_step=1e-5)
    second_gap = float(sols[1]["objective"]-sols[0]["objective"]) if len(sols) > 1 else None
    return refined, second_gap


def write_csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def main(config_path):
    cfg_path = Path(config_path).resolve()
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    out = ROOT/"outputs"/"boundary_scan"/cfg["run_id"]
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True)
    angles = np.linspace(cfg["angle_min_deg"], cfg["angle_max_deg"], cfg["angle_bins"])
    materials = {k: ROOT/v for k, v in cfg["materials"].items()}
    jc = materials["JC"]
    started = time.time()

    local = []
    for budget in cfg["total_incident_photons"]:
        for d in cfg["gold_thickness_nm"]:
            for dq in cfg["sample_delta_indices"]:
                for fs in cfg["source_reference_fractions"]:
                    local.append(local_row(cfg, angles, jc, budget, d, dq, fs))
    write_csv(out/"local_geometry.csv", local)
    print(f"local geometry rows={len(local)}", flush=True)

    t2rows = []
    for budget in cfg["total_incident_photons"]:
        for d in cfg["gold_thickness_nm"]:
            for dq in cfg["sample_delta_indices"]:
                for allocation in cfg["t2_allocations"]:
                    t2rows.append(t2_local_row(cfg, angles, jc, budget, d, dq, allocation))
    write_csv(out/"t2_local_geometry.csv", t2rows)
    print(f"T2 local geometry rows={len(t2rows)}", flush=True)

    a = cfg["anchor"]
    q = cfg["baseline_index"]+a["sample_delta_index"]
    d = a["gold_thickness_nm"]
    budget = a["total_incident_photons"]
    pair_alloc = [a["source_reference_fraction"]]
    t2_alloc = a["t2_allocation"]
    anchor_geom = local_row(cfg, angles, jc, budget, d, a["sample_delta_index"], pair_alloc[0])
    pseudo = []
    for delta in cfg["mismatch_nm"]:
        for strategy, alloc in (("T0", pair_alloc), ("T1", pair_alloc), ("T2", t2_alloc), ("O", pair_alloc)):
            y = truth_vector(cfg, angles, jc, jc, budget, alloc, q, d, delta, t2=(strategy == "T2"))
            fit, gap = deterministic_fit(strategy, y, cfg, angles, jc, budget, alloc,
                                         known_delta=delta if strategy == "O" else 0.0)
            estq = float(fit["params"][0])
            pseudo.append({
                "truth_material_source": "JC", "truth_material_target": "JC",
                "fit_material": "JC", "strategy": strategy, "true_mismatch_nm": delta,
                "true_target_angle_offset_deg": 0.0, "truth_q": q, "estimated_q": estq,
                "q_bias_riu": estq-q,
                "linear_predicted_T0_bias_riu": anchor_geom["K_delta_riu_per_nm"]*delta if strategy == "T0" else "",
                "bias_minus_linear_prediction_riu": (estq-q-anchor_geom["K_delta_riu_per_nm"]*delta) if strategy == "T0" else "",
                "objective_half_deviance": fit["objective"], "fit_iterations": fit["iterations"],
                "fit_accepted_steps": fit["accepted_steps"], "coarse_second_best_gap": gap,
                "estimated_parameters_json": json.dumps([float(x) for x in fit["params"]])
            })
        print(f"pseudo-true mismatch_nm={delta}", flush=True)
    write_csv(out/"pseudo_true_boundary.csv", pseudo)

    angle_rows = []
    for tau in cfg["target_angle_offset_deg"]:
        for strategy, alloc, aware in (("T0", pair_alloc, False), ("T1", pair_alloc, False),
                                       ("T1_angle", pair_alloc, True), ("T2", t2_alloc, False),
                                       ("T2_angle", t2_alloc, True)):
            base = "T1" if strategy.startswith("T1") else "T2" if strategy.startswith("T2") else "T0"
            y = truth_vector(cfg, angles, jc, jc, budget, alloc, q, d, 0.0, tau, t2=(base == "T2"))
            fit, gap = deterministic_fit(base, y, cfg, angles, jc, budget, alloc, fit_target_angle=aware)
            angle_rows.append({
                "strategy": strategy, "true_mismatch_nm": 0.0, "true_target_angle_offset_deg": tau,
                "truth_q": q, "estimated_q": float(fit["params"][0]),
                "q_bias_riu": float(fit["params"][0])-q,
                "linear_predicted_T0_bias_riu": anchor_geom["K_target_angle_riu_per_deg"]*tau if strategy == "T0" else "",
                "objective_half_deviance": fit["objective"], "coarse_second_best_gap": gap,
                "estimated_parameters_json": json.dumps([float(x) for x in fit["params"]])
            })
        print(f"angle-offset pseudo-true tau_deg={tau}", flush=True)
    write_csv(out/"target_angle_offset.csv", angle_rows)

    material_rows = []
    cases = []
    for name in ("Olmon_evaporated", "Rakic_LD"):
        cases.append(("public_model_mismatch", name, name))
        cases.append(("interdevice_material_change", "JC", name))
    for case, source_name, target_name in cases:
        for strategy, alloc in (("T0", pair_alloc), ("T1", pair_alloc), ("T2", t2_alloc), ("O", pair_alloc)):
            y = truth_vector(cfg, angles, materials[source_name], materials[target_name], budget, alloc,
                             q, d, 0.0, t2=(strategy == "T2"))
            fit, gap = deterministic_fit(strategy, y, cfg, angles, jc, budget, alloc, known_delta=0.0)
            material_rows.append({
                "case": case, "truth_source_material": source_name, "truth_target_material": target_name,
                "fit_material": "JC", "strategy": strategy, "true_mismatch_nm": 0.0,
                "truth_q": q, "estimated_q": float(fit["params"][0]),
                "q_bias_riu": float(fit["params"][0])-q,
                "objective_half_deviance": fit["objective"], "coarse_second_best_gap": gap,
                "estimated_parameters_json": json.dumps([float(x) for x in fit["params"]])
            })
        print(f"material case={case} source={source_name} target={target_name}", flush=True)
    write_csv(out/"material_mismatch.csv", material_rows)

    small = [r for r in pseudo if r["strategy"] == "T0" and 0 < abs(r["true_mismatch_nm"]) <= 0.2]
    pred = np.array([float(r["linear_predicted_T0_bias_riu"]) for r in small])
    obs = np.array([r["q_bias_riu"] for r in small])
    summary = {
        "status": "EXPLORATORY_COMPLETE", "run_id": cfg["run_id"], "seed": cfg["seed"],
        "elapsed_seconds": time.time()-started, "config_sha256": digest(cfg_path),
        "material_sha256": {k: digest(v) for k, v in materials.items()},
        "observation_model": {
            "interrogation": "angular response at 633 nm over 60-80 degrees",
            "angle_bins": cfg["angle_bins"],
            "photon_budget": "total across all bins and all exposures; divided by exposure allocation then evenly across bins",
            "background": "0.001 of incident photons in each bin, scaled with that exposure's photon allocation",
            "likelihood": "independent Poisson counts in the abstract bin model"
        },
        "anchor": a,
        "anchor_K_delta_riu_per_nm": anchor_geom["K_delta_riu_per_nm"],
        "anchor_K_target_angle_riu_per_deg": anchor_geom["K_target_angle_riu_per_deg"],
        "anchor_T0_vs_T1_boundary_abs_delta_nm": anchor_geom["T0_vs_T1_boundary_abs_delta_nm"],
        "small_mismatch_linear_prediction_rmse_riu": float(np.sqrt(np.mean((obs-pred)**2))),
        "small_mismatch_linear_prediction_max_abs_error_riu": float(np.max(np.abs(obs-pred))),
        "terminology_correction": "derivative norm retained is a length ratio; information fraction retained is its square",
        "interpretation_limit": "scenario-based local geometry and deterministic pseudo-true fitting; no Monte Carlo confirmation or empirical drift distribution"
    }
    (out/"summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_boundary_scan.py CONFIG")
    main(sys.argv[1])
