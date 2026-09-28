#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from calspr.optics import nk_gold_johnson, reflectance_kretschmann


def digest(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()


def cdiff(fn,x,h): return (fn(x+h)-fn(x-h))/(2*h)


def relnorm(a,b):
    den=max(float(np.linalg.norm(b)),np.finfo(float).tiny)
    return float(np.linalg.norm(a-b)/den)


def main(config_path):
    p=Path(config_path).resolve(); cfg=json.loads(p.read_text(encoding="utf-8"))
    out=ROOT/"outputs"/"derivative_validation"/cfg["run_id"]
    if out.exists(): raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True)
    angles=np.linspace(cfg["angle_min_deg"],cfg["angle_max_deg"],cfg["angle_bins"])
    rows=[]; started=time.time()
    for material,rel in cfg["materials"].items():
        source=ROOT/rel; nk=complex(nk_gold_johnson(cfg["wavelength_nm"],source))
        for d in cfg["thickness_nm"]:
            for q in cfg["analyte_index"]:
                for param,base in cfg["base_steps"].items():
                    derivs={}
                    for method in ("recursive","characteristic"):
                        fn0=lambda dd,qq,aa: reflectance_kretschmann(angles+aa,cfg["wavelength_nm"],dd,qq,source,method=method)
                        for scale in cfg["step_scales"]:
                            h=base*scale
                            if param=="refractive_index": val=cdiff(lambda x:fn0(d,x,0.0),q,h)
                            elif param=="thickness_nm": val=cdiff(lambda x:fn0(x,q,0.0),d,h)
                            else: val=cdiff(lambda x:fn0(d,q,x),0.0,h)
                            derivs[(method,scale)]=val
                    fine=cfg["step_scales"][-1]
                    rows.append({
                        "material":material,"gold_n":nk.real,"gold_k":nk.imag,
                        "thickness_nm":d,"analyte_index":q,"parameter":param,
                        "base_step":base,
                        "recursive_step1_vs_half_relative_norm":relnorm(derivs[("recursive",1.0)],derivs[("recursive",0.5)]),
                        "recursive_half_vs_quarter_relative_norm":relnorm(derivs[("recursive",0.5)],derivs[("recursive",0.25)]),
                        "solver_derivative_relative_norm_at_quarter":relnorm(derivs[("recursive",fine)],derivs[("characteristic",fine)]),
                        "recursive_derivative_norm_at_quarter":float(np.linalg.norm(derivs[("recursive",fine)])),
                        "all_finite":int(all(np.all(np.isfinite(x)) for x in derivs.values()))
                    })
    with open(out/"derivative_checks.csv","w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    summary={
        "status":"NUMERICAL_VERIFICATION_COMPLETE","run_id":cfg["run_id"],"seed":cfg["seed"],
        "elapsed_seconds":time.time()-started,"config_sha256":digest(p),
        "material_sha256":{k:digest(ROOT/v) for k,v in cfg["materials"].items()},
        "n_checks":len(rows),"all_finite":bool(all(r["all_finite"] for r in rows)),
        "max_recursive_half_vs_quarter_relative_norm":max(r["recursive_half_vs_quarter_relative_norm"] for r in rows),
        "max_solver_derivative_relative_norm_at_quarter":max(r["solver_derivative_relative_norm_at_quarter"] for r in rows),
        "scope":"Centered-derivative step convergence and agreement of two matrix-optics implementations at 633 nm across named thickness, index, and material states."
    }
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    if len(sys.argv)!=2: raise SystemExit("usage: run_derivative_validation.py CONFIG")
    main(sys.argv[1])
