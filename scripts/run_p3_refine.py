#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from calspr.inversion import fisher_scoring_refine
from run_p3_inversion import model_vector, truth_vector, target_truth, summarize


def main(config_path):
    rcfg=json.loads(Path(config_path).read_text(encoding="utf-8"))
    cfg=json.loads((ROOT/rcfg["source_config"]).read_text(encoding="utf-8"))
    source=ROOT/"outputs"/"p3_inversion"/rcfg["source_run_id"]/"replicate_results.csv"
    out=ROOT/"outputs"/"p3_refine"/rcfg["run_id"]
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True)
    with open(source,encoding="utf-8") as f:
        old=list(csv.DictReader(f))
    lookup={(r["regime"],float(r["true_transfer_mismatch_nm"]),int(r["replicate"])):r for r in old}
    gold=ROOT/cfg["material_sources"]["gold_raw"]
    angles=np.linspace(cfg["angle_min_deg"],cfg["angle_max_deg"],cfg["angle_bins"])
    rng=np.random.default_rng(cfg["seed"])
    conditions=[("absolute",0.0),("same_device_differential",0.0)]
    conditions += [(regime,m) for regime in ("transfer_naive","transfer_aware") for m in cfg["transfer_mismatch_nm"]]
    records=[]; started=time.time()
    for regime,mismatch in conditions:
        mu=truth_vector(regime,mismatch,cfg,angles,gold)
        truth=target_truth(regime,cfg)
        bounds=cfg["bounds"][regime]
        for rep in range(cfg["replicates_per_condition"]):
            y=rng.poisson(mu)
            prior=lookup[(regime,float(mismatch),rep)]
            initial=json.loads(prior["best_parameters_json"])
            refined=fisher_scoring_refine(y,lambda p:model_vector(regime,p,cfg,angles,gold),initial,bounds,
                                           max_iter=rcfg["max_iterations"],derivative_step=rcfg["unit_derivative_step"])
            target=float(refined["params"][0]); lo,hi=bounds[0]
            bfrac=cfg["fit"]["target_boundary_fraction"]
            records.append({"regime":regime,"true_transfer_mismatch_nm":mismatch,"replicate":rep,
                            "truth_target":truth,"estimated_target":target,"target_error_riu":target-truth,
                            "objective_half_deviance":refined["objective"],"optimizer_failure":int(not np.isfinite(refined["objective"])),
                            "target_boundary":int((target-lo<=bfrac*(hi-lo)) or (hi-target<=bfrac*(hi-lo))),
                            "multi_solution":int(prior["multi_solution"]),"evaluations_all_starts":int(prior["evaluations_all_starts"]),
                            "refine_iterations":refined["iterations"],"refine_accepted_steps":refined["accepted_steps"],
                            "objective_improvement":float(prior["objective_half_deviance"])-refined["objective"],
                            "best_parameters_json":json.dumps([float(x) for x in refined["params"]])})
        print(f"refined regime={regime} mismatch_nm={mismatch}",flush=True)
    with open(out/"refined_replicates.csv","w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
    grouped=[]
    for regime,mismatch in conditions:
        g=[r for r in records if r["regime"]==regime and r["true_transfer_mismatch_nm"]==mismatch]
        grouped.append({"regime":regime,"true_transfer_mismatch_nm":mismatch,**summarize(g),
                        "median_objective_improvement":float(np.median([r["objective_improvement"] for r in g]))})
    with open(out/"condition_summary.csv","w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(grouped[0])); w.writeheader(); w.writerows(grouped)
    summary={"status":"REFINEMENT_COMPLETE","run_id":rcfg["run_id"],"source_run_id":rcfg["source_run_id"],
             "elapsed_seconds":time.time()-started,"seed":cfg["seed"],"conditions":grouped,
             "purpose":"distinguish coarse optimizer artifacts from likelihood geometry; same Poisson draws regenerated from frozen seed and order"}
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    if len(sys.argv)!=2: raise SystemExit("usage: run_p3_refine.py CONFIG")
    main(sys.argv[1])

