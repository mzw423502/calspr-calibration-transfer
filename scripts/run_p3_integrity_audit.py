#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from run_p3_inversion import summarize, truth_vector


def read_csv(path):
    with open(path,encoding="utf-8") as f: return list(csv.DictReader(f))


def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()


def main(config_path):
    ap=Path(config_path).resolve(); acfg=json.loads(ap.read_text(encoding="utf-8"))
    cfg=json.loads((ROOT/acfg["source_config"]).read_text(encoding="utf-8"))
    out=ROOT/"outputs"/"integrity_audit"/acfg["run_id"]
    if out.exists(): raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True)
    coarse=read_csv(ROOT/acfg["coarse_records"]); refined=read_csv(ROOT/acfg["refined_records"])
    published=read_csv(ROOT/acfg["refined_summary"])
    key=lambda r:(r["regime"],float(r["true_transfer_mismatch_nm"]),int(r["replicate"]))
    ck={key(r) for r in coarse}; rk={key(r) for r in refined}
    conditions=[("absolute",0.0),("same_device_differential",0.0)]
    conditions += [(regime,float(m)) for regime in ("transfer_naive","transfer_aware") for m in cfg["transfer_mismatch_nm"]]
    recomputed=[]
    for regime,mismatch in conditions:
        g=[]
        for r in refined:
            if r["regime"]==regime and float(r["true_transfer_mismatch_nm"])==mismatch:
                z=dict(r)
                for fld in ("target_error_riu","optimizer_failure","target_boundary","multi_solution","evaluations_all_starts"):
                    z[fld]=float(z[fld])
                g.append(z)
        recomputed.append({"regime":regime,"true_transfer_mismatch_nm":mismatch,**summarize(g)})
    pub={(r["regime"],float(r["true_transfer_mismatch_nm"])):r for r in published}
    fields=["bias_riu","rmse_riu","mae_riu","p95_absolute_error_riu","optimizer_failure_rate","target_boundary_rate","multi_solution_rate","median_evaluations_all_starts"]
    diffs=[]
    for r in recomputed:
        p=pub[(r["regime"],r["true_transfer_mismatch_nm"])]
        for f in fields: diffs.append(abs(float(r[f])-float(p[f])))

    angles=np.linspace(cfg["angle_min_deg"],cfg["angle_max_deg"],cfg["angle_bins"])
    gold=ROOT/acfg["material"]
    rng=np.random.default_rng(cfg["seed"])
    fingerprints=[]
    for regime,mismatch in conditions:
        mu=truth_vector(regime,mismatch,cfg,angles,gold)
        h=hashlib.sha256()
        for rep in range(cfg["replicates_per_condition"]):
            y=rng.poisson(mu); h.update(np.asarray(y,dtype=np.int64).tobytes(order="C"))
        fingerprints.append({"regime":regime,"true_transfer_mismatch_nm":mismatch,
                             "replicates":cfg["replicates_per_condition"],"poisson_counts_sha256":h.hexdigest(),
                             "observation_bins_per_replicate":len(mu)})
    with open(out/"recomputed_condition_summary.csv","w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(recomputed[0])); w.writeheader(); w.writerows(recomputed)
    with open(out/"poisson_reconstruction_fingerprints.csv","w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(fingerprints[0])); w.writeheader(); w.writerows(fingerprints)
    condition_counts={f"{a}|{b}":sum(1 for r in refined if r["regime"]==a and float(r["true_transfer_mismatch_nm"])==b) for a,b in conditions}
    report={
        "status":"PASS" if ck==rk and len(refined)==240 and max(diffs)<=acfg["numeric_tolerance"] else "FAIL",
        "run_id":acfg["run_id"],"numpy_version":np.__version__,"rng":"numpy.random.default_rng (PCG64)",
        "source_config_sha256":sha(ROOT/acfg["source_config"]),"coarse_records_sha256":sha(ROOT/acfg["coarse_records"]),
        "refined_records_sha256":sha(ROOT/acfg["refined_records"]),"published_summary_sha256":sha(ROOT/acfg["refined_summary"]),
        "condition_count":len(conditions),"refined_record_count":len(refined),"coarse_record_count":len(coarse),
        "condition_counts":condition_counts,"coarse_refined_keys_identical":ck==rk,
        "max_absolute_summary_difference":max(diffs),"numeric_tolerance":acfg["numeric_tolerance"],
        "summary_definitions":{"bias":"mean target error in RIU","rmse":"square root of mean squared target error in RIU","mae":"mean absolute target error in RIU","p95_absolute_error":"numpy linear-interpolated 0.95 quantile of absolute target error in RIU"},
        "reconstruction_limit":"Fingerprints prove the counts reconstructed by this recorded generator, version, seed, condition order, and call order. Counts were not stored by the original runs, so this audit cannot independently hash the transient arrays used then; source-code and key-order agreement are additional evidence, not a stored-count identity certificate."
    }
    (out/"report.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))
    if report["status"]!="PASS": raise SystemExit(2)


if __name__=="__main__":
    if len(sys.argv)!=2: raise SystemExit("usage: run_p3_integrity_audit.py CONFIG")
    main(sys.argv[1])
