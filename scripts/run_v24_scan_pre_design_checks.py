#!/usr/bin/env python3
"""CalSPR V2.4 bounded information-source and nested-threshold scan."""
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
from calspr.v24_design import (
    ResponseScenario,
    angle_centres,
    bare_prism_reflectance,
    finite_difference_columns,
    gold_reflectance,
    information_diagnostics,
    physical_windows,
    poisson_mean,
    response_scenarios,
)


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"no rows for {path}")
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def scenario_by_name(cfg, name):
    return next(s for s in response_scenarios(cfg) if s.name == name)


def _gold_exposure(cfg, centres, gold, scenario, photons, d, q, a0, lg):
    r = gold_reflectance(centres + a0, cfg["wavelength_nm"], d, q, gold, scenario)
    return poisson_mean(r, photons, cfg["background_fraction"], lg)


def _prism_exposure(cfg, centres, scenario, photons, q, a0, lg):
    r = bare_prism_reflectance(centres + a0, cfg["wavelength_nm"], q, scenario)
    return poisson_mean(r, photons, cfg["background_fraction"], lg)


def angular_segments(strategy, windows, critical_weight=0.5):
    if strategy == "CS":
        return [(windows["C"], float(critical_weight)),
                (windows["S"], 1.0 - float(critical_weight))]
    if strategy == "P":
        return [(windows["F"], 1.0)]
    return [(windows[strategy], 1.0)]


def design_model(cfg, centres, gold, scenario, windows, strategy, ref_fraction,
                 total_photons, critical_weight=0.5):
    segments = angular_segments(strategy, windows, critical_weight)
    q0 = float(cfg["baseline_index"])
    if strategy == "P":
        def fn(p):
            q, a0, lg = p
            out = []
            for mask, aw in segments:
                if aw <= 0:
                    continue
                ang = centres[mask]
                out.append(_prism_exposure(cfg, ang, scenario, total_photons*ref_fraction*aw,
                                           q0, a0, lg))
                out.append(_prism_exposure(cfg, ang, scenario, total_photons*(1-ref_fraction)*aw,
                                           q, a0, lg))
            return np.concatenate(out)
        p = np.array([q0+cfg["sample_delta_index"], 0.0, 0.0])
        h = cfg["finite_difference_steps"]
        steps = [h["refractive_index"], h["angle_zero_deg"], h["log_gain"]]
        labels = ["q", "angle_zero_deg", "log_gain"]
    else:
        def fn(p):
            q, d, a0, lg = p
            out = []
            for mask, aw in segments:
                if aw <= 0:
                    continue
                ang = centres[mask]
                out.append(_gold_exposure(cfg, ang, gold, scenario,
                                          total_photons*ref_fraction*aw, d, q0, a0, lg))
                out.append(_gold_exposure(cfg, ang, gold, scenario,
                                          total_photons*(1-ref_fraction)*aw, d, q, a0, lg))
            return np.concatenate(out)
        p = np.array([q0+cfg["sample_delta_index"], cfg["nominal_thickness_nm"], 0.0, 0.0])
        h = cfg["finite_difference_steps"]
        steps = [h["refractive_index"], h["thickness_nm"], h["angle_zero_deg"], h["log_gain"]]
        labels = ["q", "thickness_nm", "angle_zero_deg", "log_gain"]
    return fn, p, steps, labels


def attribution_rows(cfg, centres, gold, scenario, windows):
    rows = []
    nfull = len(centres)
    budget = float(cfg["equal_budget_total_incident_photons"])
    ref_fraction = 0.5
    for design in ("F", "no_C", "no_S", "C", "S", "CS", "O"):
        mask = windows[design]
        used_fraction = float(np.sum(mask)) / nfull
        uniform_cs_weight = (float(np.sum(windows["C"])) / float(np.sum(mask))
                             if design == "CS" else 0.5)
        # Discarded bins are not reallocated: preserve full-design photons/bin.
        fn, p, steps, labels = design_model(
            cfg, centres, gold, scenario, windows, design, ref_fraction,
            budget * used_fraction, uniform_cs_weight,
        )
        _, cols = finite_difference_columns(fn, p, steps)
        info = information_diagnostics(cols)
        rows.append({
            "design": design,
            "angular_response": scenario.name,
            "bin_count": int(np.sum(mask)),
            "used_photon_fraction": used_fraction,
            "total_incident_photons_used": budget * used_fraction,
            "discarded_photons_reallocated": 0,
            "parameter_labels": json.dumps(labels),
            "rank": info["rank"],
            "parameter_count": info["parameter_count"],
            "condition_number_whitened_jacobian": info["condition_number_whitened_jacobian"],
            "effective_information": info["effective_information"],
            "local_sd_riu": info["local_sd"],
            "identifiable": int(info["rank"] == info["parameter_count"]),
        })

    # Full-design projected-score contributions. These sum to the jointly
    # profiled F information and are not separately re-profiled region infos.
    fn, p, steps, _ = design_model(cfg, centres, gold, scenario, windows, "F",
                                   ref_fraction, budget)
    _, cols = finite_difference_columns(fn, p, steps)
    info = information_diagnostics(cols)
    rq2 = info["projected_target_score"] ** 2
    for name in ("C", "S", "O"):
        mask = windows[name]
        # model order is all reference rows followed by all sample rows
        contribution = float(np.sum(rq2[:nfull][mask]) + np.sum(rq2[nfull:][mask]))
        rows.append({
            "design": f"F_projected_score_contribution_{name}",
            "angular_response": scenario.name,
            "bin_count": int(np.sum(mask)),
            "used_photon_fraction": 1.0,
            "total_incident_photons_used": budget,
            "discarded_photons_reallocated": 0,
            "parameter_labels": "joint F projection; contribution only",
            "rank": info["rank"],
            "parameter_count": info["parameter_count"],
            "condition_number_whitened_jacobian": info["condition_number_whitened_jacobian"],
            "effective_information": contribution,
            "local_sd_riu": "",
            "identifiable": 1,
        })
    return rows


def equal_budget_rows(cfg, centres, gold, scenario, windows):
    rows = []
    budget = float(cfg["equal_budget_total_incident_photons"])
    for strategy in ("F", "C", "S", "CS", "P"):
        cweights = cfg["cs_critical_weights"] if strategy == "CS" else [""]
        for ref_fraction in cfg["reference_fractions"]:
            for cw in cweights:
                critical_weight = 0.5 if cw == "" else float(cw)
                fn, p, steps, labels = design_model(
                    cfg, centres, gold, scenario, windows, strategy,
                    float(ref_fraction), budget, critical_weight,
                )
                _, cols = finite_difference_columns(fn, p, steps)
                info = information_diagnostics(cols)
                rows.append({
                    "strategy": strategy,
                    "angular_response": scenario.name,
                    "total_incident_photons": budget,
                    "reference_fraction": ref_fraction,
                    "sample_fraction": 1.0-float(ref_fraction),
                    "cs_critical_weight": cw,
                    "cs_spr_weight": "" if cw == "" else 1.0-float(cw),
                    "bin_count_C": int(np.sum(windows["C"])) if strategy == "CS" else "",
                    "bin_count_S": int(np.sum(windows["S"])) if strategy == "CS" else "",
                    "bin_count_total": len(centres) if strategy == "P" else int(np.sum(windows[strategy])),
                    "parameter_labels": json.dumps(labels),
                    "rank": info["rank"],
                    "parameter_count": info["parameter_count"],
                    "condition_number_whitened_jacobian": info["condition_number_whitened_jacobian"],
                    "effective_information": info["effective_information"],
                    "local_sd_riu": info["local_sd"],
                    "identifiable": int(info["rank"] == info["parameter_count"]),
                    "selected_within_strategy": 0,
                    "threshold_applicable": "NO" if strategy == "P" else "not_evaluated_here",
                })
    for strategy in ("F", "C", "S", "CS", "P"):
        candidates = [r for r in rows if r["strategy"] == strategy and r["identifiable"]]
        if candidates:
            min(candidates, key=lambda r: float(r["local_sd_riu"]))["selected_within_strategy"] = 1
    return rows


def transfer_full_model(cfg, centres, gold, scenario, budget, thickness):
    q0 = float(cfg["baseline_index"])

    def fn(p):
        q, d, a0, lg, delta = p
        return np.concatenate([
            _gold_exposure(cfg, centres, gold, scenario, budget*0.5, d, q0, a0, lg),
            _gold_exposure(cfg, centres, gold, scenario, budget*0.5, d+delta, q, a0, lg),
        ])

    p = np.array([q0+cfg["sample_delta_index"], thickness, 0.0, 0.0, 0.0])
    h = cfg["finite_difference_steps"]
    steps = [h["refractive_index"], h["thickness_nm"], h["angle_zero_deg"],
             h["log_gain"], h["thickness_nm"]]
    return fn, p, steps


def threshold_row(cfg, centres, gold, scenario, budget, thickness, step_multiplier=1.0,
                  phase_fraction=0.0):
    fn, p, steps = transfer_full_model(cfg, centres, gold, scenario, budget, thickness)
    steps = [float(h)*float(step_multiplier) for h in steps]
    _, cols = finite_difference_columns(fn, p, steps)
    fisher = cols.T @ cols
    a = fisher[:4, :4]
    b = fisher[:4, 4]
    c = float(fisher[4, 4])
    a_inv = np.linalg.inv(a)
    full_inv = np.linalg.inv(fisher)
    k = float((a_inv @ b)[0])
    v0 = float(a_inv[0, 0])
    v1 = float(full_inv[0, 0])
    schur = float(c - b @ a_inv @ b)
    rhs = float(k*k/schur)
    dstar_general = math.sqrt(max(0.0, v1-v0))/abs(k) if k != 0 else math.nan
    dstar_schur = 1.0/math.sqrt(schur) if schur > 0 else math.nan
    norms = np.linalg.norm(cols, axis=0)
    scaled_cols = cols / norms
    scaled_s = np.linalg.svd(scaled_cols, compute_uv=False)
    rank = int(np.sum(scaled_s > scaled_s[0]*1e-10))
    return {
        "anchor_id": f"d{thickness:g}_N{budget:.0e}",
        "thickness_nm": thickness,
        "total_incident_photons": budget,
        "angular_response": scenario.name,
        "phase_fraction_of_spacing": phase_fraction,
        "finite_difference_step_multiplier": step_multiplier,
        "full_fisher_rank": rank,
        "full_parameter_count": 5,
        "A_condition_number": float(np.linalg.cond(a)),
        "full_fisher_condition_number": float(np.linalg.cond(fisher)),
        "scaled_whitened_jacobian_condition_number": float(scaled_s[0]/scaled_s[-1]),
        "K_delta_riu_per_nm": k,
        "v0_riu2": v0,
        "v1_riu2": v1,
        "v1_minus_v0_riu2": v1-v0,
        "S_delta_per_nm2": schur,
        "K2_over_S_riu2": rhs,
        "variance_identity_relative_error": abs((v1-v0)-rhs)/max(abs(v1-v0), abs(rhs), 1e-300),
        "delta_star_general_nm": dstar_general,
        "delta_star_schur_nm": dstar_schur,
        "threshold_identity_relative_error": abs(dstar_general-dstar_schur)/max(abs(dstar_schur), 1e-300),
        "regular_interior_identifiable": int(rank == 5 and schur > 0),
    }


def numerical_checks(cfg, gold, scenarios):
    rows = []
    main = scenario_by_name(cfg, cfg["main_angular_response"])
    for scenario in scenarios:
        for phase in cfg["grid_phase_fractions"]:
            centres = angle_centres(cfg, phase)
            row = threshold_row(cfg, centres, gold, scenario, 1e7, 50.0,
                                step_multiplier=1.0, phase_fraction=phase)
            row["check_type"] = "response_and_grid_phase"
            rows.append(row)
    for mult in cfg["finite_difference_step_multipliers"]:
        centres = angle_centres(cfg, 0.0)
        row = threshold_row(cfg, centres, gold, main, 1e7, 50.0,
                            step_multiplier=mult, phase_fraction=0.0)
        row["check_type"] = "finite_difference_step"
        rows.append(row)
    for order in cfg["quadrature_orders_to_check"]:
        scenario = ResponseScenario(main.name+f"_quad{order}", main.kind, main.width_deg, int(order))
        centres = angle_centres(cfg, 0.0)
        row = threshold_row(cfg, centres, gold, scenario, 1e7, 50.0)
        row["check_type"] = "quadrature_order"
        rows.append(row)
    return rows


def main(config_path):
    started = time.time()
    cfg_path = Path(config_path).resolve()
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    out = ROOT / "outputs" / "v24" / cfg["run_id"]
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True)
    gold = ROOT / cfg["material"]
    centres = angle_centres(cfg)
    spacing = float(centres[1]-centres[0])
    scenarios = response_scenarios(cfg)
    main_scenario = scenario_by_name(cfg, cfg["main_angular_response"])
    windows = physical_windows(cfg, centres, gold)

    window_rows = []
    for name in ("C", "S", "O", "CS", "F", "no_C", "no_S"):
        ang = centres[windows[name]]
        window_rows.append({
            "window": name, "bin_count": len(ang), "angle_min_deg": float(np.min(ang)),
            "angle_max_deg": float(np.max(ang)), "definition_basis": cfg["window_definition"]["basis"],
            "critical_angle_deg": windows["critical_angle_deg"],
            "spr_minimum_deg": windows["spr_minimum_deg"],
        })
    write_csv(out / "V24_WINDOWS.csv", window_rows)

    attribution = attribution_rows(cfg, centres, gold, main_scenario, windows)
    write_csv(out / "V24_INFORMATION_ATTRIBUTION.csv", attribution)
    equal = equal_budget_rows(cfg, centres, gold, main_scenario, windows)
    write_csv(out / "V24_EQUAL_BUDGET_COMPARISON.csv", equal)

    threshold = []
    for d in cfg["anchor_thickness_nm"]:
        for budget in cfg["anchor_total_incident_photons"]:
            threshold.append(threshold_row(cfg, centres, gold, main_scenario, budget, d))
    write_csv(out / "V24_THRESHOLD_PREDICTIONS.csv", threshold)
    checks = numerical_checks(cfg, gold, scenarios)
    write_csv(out / "V24_NUMERICAL_CHECKS.csv", checks)

    # Exact legacy point-response anchor re-run, compared with the recorded V2.2
    # values only as a regression check (rounded manuscript values are labelled).
    legacy = scenario_by_name(cfg, "legacy_point")
    legacy_row = threshold_row(cfg, centres, gold, legacy, 1e7, 50.0)
    expected = cfg["legacy_anchor_expected"]
    legacy_comparison = {
        "computed": legacy_row,
        "recorded_v22_rounded_or_exact": expected,
        "relative_difference": {
            "K_delta_riu_per_nm": abs(legacy_row["K_delta_riu_per_nm"]-expected["K_delta_riu_per_nm"])/abs(expected["K_delta_riu_per_nm"]),
            "delta_star_nm": abs(legacy_row["delta_star_schur_nm"]-expected["delta_star_nm"])/abs(expected["delta_star_nm"]),
        },
        "note": "V2.2 information values in the manuscript are rounded; delta_star is the recorded exact value.",
    }
    (out / "V24_LEGACY_POINT_REPRODUCTION.json").write_text(
        json.dumps(legacy_comparison, indent=2)+"\n", encoding="utf-8")

    selected = {r["strategy"]: r for r in equal if r["selected_within_strategy"] == 1}
    summary = {
        "status": "V24_BOUNDED_SCAN_COMPLETE",
        "run_id": cfg["run_id"],
        "seed": cfg["seed"],
        "elapsed_seconds": time.time()-started,
        "config_sha256": digest(cfg_path),
        "material_sha256": digest(gold),
        "angle_spacing_deg": spacing,
        "point_semantics": "legacy centre-point reflectance; independent Poisson means",
        "finite_response_semantics": "normalized top-hat angular integration; named numerical scenarios without empirical instrument attribution",
        "window_definition": {k: windows[k] for k in ("critical_angle_deg", "spr_minimum_deg")},
        "selected_equal_budget_designs": selected,
        "threshold_anchor_count": len(threshold),
        "all_threshold_identities_pass_1e-6": all(r["variance_identity_relative_error"] < 1e-6 and r["threshold_identity_relative_error"] < 1e-6 for r in threshold),
        "interpretation_limit": cfg["interpretation_limit"],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_v24_scan.py CONFIG")
    main(sys.argv[1])
