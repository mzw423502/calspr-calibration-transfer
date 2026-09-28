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
from calspr.optics import nk_gold_johnson, reflectance_kretschmann


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main(config_path):
    cfg_path = Path(config_path).resolve()
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    out = ROOT / "outputs" / "p4_materials" / cfg["run_id"]
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True)
    angles = np.linspace(cfg["angle_min_deg"], cfg["angle_max_deg"], cfg["angle_bins"])
    info, resonance, hashes = [], [], {}
    start = time.time()
    for material, relpath in cfg["materials"].items():
        source = ROOT / relpath
        hashes[material] = digest(source)
        nk = complex(nk_gold_johnson(cfg["wavelength_nm"], source))
        for d in cfg["gold_thickness_nm"]:
            for state, q in (("baseline", cfg["baseline_index"]), ("sample", cfg["baseline_index"]+cfg["sample_delta_index"])):
                r = reflectance_kretschmann(angles, cfg["wavelength_nm"], d, q, source)
                j = int(np.argmin(r))
                resonance.append({"material":material,"gold_n":nk.real,"gold_k":nk.imag,"thickness_nm":d,"state":state,
                                  "resonance_angle_deg":float(angles[j]),"minimum_reflectance":float(r[j])})
            for regime in ("absolute","same_device_differential","cross_device_transfer"):
                deriv, labels, _, _ = scenario_derivatives(cfg, d, regime, source)
                raw, _, sep = effective_information(deriv)
                hybrid, _, hsep = effective_information(deriv, prior_rows=calibration_prior_rows(cfg, labels))
                info.append({"material":material,"gold_n":nk.real,"gold_k":nk.imag,"thickness_nm":d,"regime":regime,
                             "observation_only_local_bound_riu":raw**-0.5 if raw>0 else None,
                             "observation_separation_fraction":sep,
                             "scenario_prior_local_bound_riu":hybrid**-0.5 if hybrid>0 else None,
                             "scenario_prior_separation_fraction":hsep})
    with open(out/"material_information.csv","w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(info[0])); w.writeheader(); w.writerows(info)
    with open(out/"material_resonance.csv","w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(resonance[0])); w.writeheader(); w.writerows(resonance)
    best=[]
    for material in cfg["materials"]:
        for regime in ("absolute","same_device_differential","cross_device_transfer"):
            rows=[x for x in info if x["material"]==material and x["regime"]==regime]
            row=min(rows,key=lambda x:x["observation_only_local_bound_riu"])
            best.append({"material":material,"regime":regime,"best_thickness_nm":row["thickness_nm"],
                         "best_observation_only_local_bound_riu":row["observation_only_local_bound_riu"],
                         "separation_fraction_at_best":row["observation_separation_fraction"]})
    summary={"status":"EXPLORATORY_COMPLETE","run_id":cfg["run_id"],"elapsed_seconds":time.time()-start,
             "seed":cfg["seed"],"config_sha256":digest(cfg_path),"material_sha256":hashes,"best_by_material_and_regime":best,
             "interpretation_limit":"differences across traceable optical-constant sources are a sensitivity analysis, not a probability distribution or external validation"}
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    if len(sys.argv)!=2:
        raise SystemExit("usage: run_p4_materials.py CONFIG")
    main(sys.argv[1])

