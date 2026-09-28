#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from calspr.information import calibration_prior_rows, effective_information, scenario_derivatives
from calspr.optics import nk_bk7, nk_gold_johnson, reflectance_kretschmann


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main(config_path):
    cfg_path = Path(config_path).resolve()
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    gold = ROOT / cfg["material_sources"]["gold_raw"]
    angles = np.linspace(cfg["angle_min_deg"], cfg["angle_max_deg"], cfg["angle_bins"])
    rows = []
    started = time.time()
    for d in cfg["gold_thickness_nm"]:
        for q_label, q in (("baseline", cfg["baseline_index"]), ("sample", cfg["baseline_index"]+cfg["sample_delta_index"])):
            r = reflectance_kretschmann(angles, cfg["wavelength_nm"], d, q, gold)
            idx = int(np.argmin(r))
            rows.append({"thickness_nm": d, "state": q_label, "analyte_index": q,
                         "resonance_angle_deg": float(angles[idx]), "minimum_reflectance": float(r[idx])})

    info_rows = []
    regimes = ("absolute", "same_device_differential", "cross_device_transfer")
    for d in cfg["gold_thickness_nm"]:
        for regime in regimes:
            deriv, labels, _, _ = scenario_derivatives(cfg, d, regime, gold)
            raw_info, raw_rank, raw_sep = effective_information(deriv)
            priors = calibration_prior_rows(cfg, labels)
            hybrid_info, hybrid_rank, hybrid_sep = effective_information(deriv, prior_rows=priors)
            info_rows.append({
                "thickness_nm": d, "regime": regime,
                "target": labels[0], "nuisance_count": len(labels)-1,
                "observation_only_effective_information": raw_info,
                "observation_only_local_bound_riu": (raw_info**-0.5 if raw_info > 0 else None),
                "observation_nuisance_rank": raw_rank,
                "observation_separation_fraction": raw_sep,
                "scenario_prior_effective_information": hybrid_info,
                "scenario_prior_local_bound_riu": (hybrid_info**-0.5 if hybrid_info > 0 else None),
                "scenario_prior_nuisance_rank": hybrid_rank,
                "scenario_prior_separation_fraction": hybrid_sep
            })

    out = ROOT / "outputs" / "p0_p2"
    out.mkdir(parents=True, exist_ok=True)
    with open(out/"resonance_scan.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with open(out/"information_scan.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(info_rows[0])); w.writeheader(); w.writerows(info_rows)

    feature_rank = 1
    summary = {
        "status": "EXPLORATORY_COMPLETE",
        "elapsed_seconds": time.time()-started,
        "seed": cfg["seed"],
        "config_sha256": sha256(cfg_path),
        "gold_source_sha256": sha256(gold),
        "materials": {"bk7_n_at_wavelength": float(nk_bk7(cfg["wavelength_nm"]).real),
                      "gold_nk_at_wavelength": str(complex(nk_gold_johnson(cfg["wavelength_nm"], gold)))},
        "single_peak_identifiability": {"observations": 1, "unknowns_example": 3,
                                        "maximum_jacobian_rank": feature_rank,
                                        "conclusion": "one peak location cannot locally identify q, thickness, and angle zero"},
        "files": ["resonance_scan.csv", "information_scan.csv"]
    }
    (out/"summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_p0_p2.py CONFIG")
    main(sys.argv[1])

