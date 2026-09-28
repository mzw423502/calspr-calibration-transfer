#!/usr/bin/env python3
"""Freeze V2.4 confirmation anchors from the completed derivative-only scan."""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main(scan_run_id, output_config):
    scan_dir = ROOT / "outputs" / "v24" / scan_run_id
    predictions = scan_dir / "V24_THRESHOLD_PREDICTIONS.csv"
    summary = scan_dir / "summary.json"
    if not predictions.exists() or not summary.exists():
        raise FileNotFoundError("completed V2.4 scan outputs are required before freezing")
    scan_state = json.loads(summary.read_text(encoding="utf-8"))
    if scan_state.get("status") != "V24_BOUNDED_SCAN_COMPLETE":
        raise RuntimeError("scan status is not complete")
    anchors = []
    with open(predictions, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if int(row["regular_interior_identifiable"]) != 1:
                raise RuntimeError(f"unidentifiable anchor: {row['anchor_id']}")
            anchors.append({
                "anchor_id": row["anchor_id"],
                "thickness_nm": float(row["thickness_nm"]),
                "total_incident_photons": float(row["total_incident_photons"]),
                "K_delta_riu_per_nm": float(row["K_delta_riu_per_nm"]),
                "v0_riu2": float(row["v0_riu2"]),
                "v1_riu2": float(row["v1_riu2"]),
                "S_delta_per_nm2": float(row["S_delta_per_nm2"]),
                "delta_pred_nm": float(row["delta_star_schur_nm"]),
            })
    if len(anchors) != 9:
        raise RuntimeError(f"expected 9 anchors, found {len(anchors)}")
    cfg = {
        "run_id": "v24_confirmation_20260914_c1",
        "run_class": "frozen_v24_paired_confirmation",
        "selection_source": str(predictions.relative_to(ROOT)),
        "selection_source_sha256": digest(predictions),
        "scan_summary_sha256": digest(summary),
        "selection_rule": "Nine precomputed anchors crossed with rho in {0, +/-0.5, +/-1, +/-2}; no noisy result entered prediction or condition selection.",
        "seed": 2026091402,
        "replicates_per_condition": 100,
        "worker_processes": 4,
        "wavelength_nm": 633.0,
        "angle_min_deg": 60.0,
        "angle_max_deg": 80.0,
        "angle_bins": 201,
        "baseline_index": 1.3330,
        "sample_delta_index": 0.0050,
        "background_fraction": 0.001,
        "angular_response": {"name": "top_hat_0p10deg_main", "kind": "top_hat", "width_deg": 0.10, "quadrature_order": 16},
        "relative_magnitudes": [0.0, -0.5, 0.5, -1.0, 1.0, -2.0, 2.0],
        "anchors": anchors,
        "bounds": {
            "q": [1.325, 1.350],
            "source_thickness_nm": [35.0, 65.0],
            "mismatch_nm": [-6.0, 6.0],
            "angle_zero_deg": [-0.10, 0.10],
            "log_gain": [-0.10, 0.10]
        },
        "fit": {
            "halton_starts": 4,
            "initial_step_fraction": 0.20,
            "step_shrink": 0.5,
            "max_levels": 9,
            "passes_per_level": 2,
            "refine_max_iterations": 60,
            "refine_unit_derivative_step": 0.00001,
            "multi_solution_deviance_delta": 2.0,
            "multi_solution_target_separation_riu": 0.0005,
            "target_boundary_fraction": 0.001
        },
        "extension_rule": "Do not extend by default. Extend the entire 9-anchor group from 100 to 200 replicates only if any |rho|=0.5 or 2 paired interval crosses zero AND its MCSE exceeds 0.25 times the absolute local predicted D; never extend to resolve rho=1, where the local predicted difference is zero.",
        "material": "data/raw/materials/Au_Johnson.yml",
        "interpretation_limit": "Paired internal Monte Carlo under a named finite-angular-response scenario; not empirical validation or a measured fabrication threshold."
    }
    target = Path(output_config)
    if not target.is_absolute():
        target = ROOT / target
    if target.exists():
        raise FileExistsError(f"refusing to overwrite {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(cfg, indent=2)+"\n", encoding="utf-8")
    benchmark = json.loads(json.dumps(cfg))
    benchmark["run_id"] = "v24_confirmation_benchmark_20260914_b1"
    benchmark["run_class"] = "v24_confirmation_runtime_memory_benchmark"
    benchmark["replicates_per_condition"] = 2
    benchmark["worker_processes"] = 1
    benchmark["relative_magnitudes"] = [0.0]
    benchmark["anchors"] = [next(a for a in anchors if a["anchor_id"] == "d50_N1e+07")]
    benchmark["interpretation_limit"] = "Runtime and memory benchmark only; not inferential evidence."
    benchmark_target = target.with_name("v24_confirmation_benchmark.json")
    if benchmark_target.exists():
        raise FileExistsError(f"refusing to overwrite {benchmark_target}")
    benchmark_target.write_text(json.dumps(benchmark, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"status": "FROZEN", "config": str(target), "sha256": digest(target),
                      "anchor_count": len(anchors), "benchmark_config": str(benchmark_target),
                      "benchmark_sha256": digest(benchmark_target)}, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: freeze_v24_confirmation.py SCAN_RUN_ID OUTPUT_CONFIG")
    main(sys.argv[1], sys.argv[2])
