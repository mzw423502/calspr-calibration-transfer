#!/usr/bin/env python3
"""Bounded noisy confirmation of V2.4 equal-budget angular designs."""
from __future__ import annotations
import csv, hashlib, json, math, sys, time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"scripts"))
from calspr.inversion import bounded_pattern_fit, fisher_scoring_refine, poisson_deviance_half
from calspr.v24_design import (ResponseScenario, angle_centres, bare_prism_reflectance,
    finite_difference_columns, gold_reflectance, information_diagnostics, physical_windows,
    poisson_mean)

def digest(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for block in iter(lambda:f.read(1<<20),b""): h.update(block)
    return h.hexdigest()

def build_model(cfg, state, strategy):
    centres=angle_centres(cfg); gold=ROOT/cfg["material"]
    a=cfg["angular_response"]; scenario=ResponseScenario(a["name"],a["kind"],a["width_deg"],a["quadrature_order"])
    windows=physical_windows(cfg,centres,gold); mask=windows["F" if strategy=="P" else strategy]
    ang=centres[mask]; fref=cfg["reference_fraction"][strategy]; total=cfg["total_incident_photons"]
    q0=cfg["baseline_index"]
    if strategy=="P":
        def fn(p):
            q,a0,lg=p
            rr=bare_prism_reflectance(ang+a0,cfg["wavelength_nm"],q0,scenario)
            rs=bare_prism_reflectance(ang+a0,cfg["wavelength_nm"],q,scenario)
            return np.concatenate([poisson_mean(rr,total*fref,cfg["background_fraction"],lg),poisson_mean(rs,total*(1-fref),cfg["background_fraction"],lg)])
        p=np.array([q0+state["sample_delta_index"],0.0,0.0]); b=cfg["bounds"]
        bounds=[b["q"],b["angle_zero_deg"],b["log_gain"]]
        h=cfg["finite_difference_steps"]; steps=[h["refractive_index"],h["angle_zero_deg"],h["log_gain"]]
    else:
        def fn(p):
            q,d,a0,lg=p
            rr=gold_reflectance(ang+a0,cfg["wavelength_nm"],d,q0,gold,scenario)
            rs=gold_reflectance(ang+a0,cfg["wavelength_nm"],d,q,gold,scenario)
            return np.concatenate([poisson_mean(rr,total*fref,cfg["background_fraction"],lg),poisson_mean(rs,total*(1-fref),cfg["background_fraction"],lg)])
        p=np.array([q0+state["sample_delta_index"],state["thickness_nm"],0.0,0.0]); b=cfg["bounds"]
        bounds=[b["q"],b["thickness_nm"],b["angle_zero_deg"],b["log_gain"]]
        h=cfg["finite_difference_steps"]; steps=[h["refractive_index"],h["thickness_nm"],h["angle_zero_deg"],h["log_gain"]]
    return fn,p,bounds,steps,len(ang)

def task(args):
    cfg,state,strategy,rep=args; fn,ptrue,bounds,steps,bins=build_model(cfg,state,strategy)
    sid=cfg["states"].index(state); mid=cfg["strategies"].index(strategy)
    stream=[cfg["seed"],sid,mid,rep]; rng=np.random.default_rng(np.random.SeedSequence(stream)); y=rng.poisson(fn(ptrue))
    sols=bounded_pattern_fit(lambda p:poisson_deviance_half(y,fn(p)),bounds,cfg["fit"]); best=sols[0]
    fit=fisher_scoring_refine(y,fn,best["params"],bounds,max_iter=cfg["fit"]["refine_max_iterations"],derivative_step=cfg["fit"]["refine_unit_derivative_step"])
    qhat=float(fit["params"][0]); qtrue=float(ptrue[0]); qfrac=(qhat-bounds[0][0])/(bounds[0][1]-bounds[0][0])
    multi=any(s["objective"]-best["objective"]<=cfg["fit"]["multi_solution_deviance_delta"] and abs(float(s["params"][0])-float(best["params"][0]))>=cfg["fit"]["multi_solution_target_separation_riu"] for s in sols[1:])
    err=qhat-qtrue
    return {"state_id":state["state_id"],"strategy":strategy,"replicate":rep,"noise_stream_json":json.dumps(stream),"thickness_nm":state["thickness_nm"],"sample_delta_index":state["sample_delta_index"],"angle_bin_count":bins,"total_incident_photons":cfg["total_incident_photons"],"reference_fraction":cfg["reference_fraction"][strategy],"true_q":qtrue,"estimated_q":qhat,"error_riu":err,"squared_error_riu2":err*err,"objective_half_deviance":float(fit["objective"]),"optimizer_failure":int(not np.isfinite(fit["objective"])),"target_boundary":int(qfrac<=cfg["fit"]["target_boundary_fraction"] or qfrac>=1-cfg["fit"]["target_boundary_fraction"]),"multi_solution":int(multi),"estimated_parameters_json":json.dumps([float(x) for x in fit["params"]])}

def write_csv(path,rows):
    keys=[]
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with open(path,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=keys); w.writeheader(); w.writerows(rows)

def main(config_path):
    started=time.time(); cp=Path(config_path).resolve(); cfg=json.loads(cp.read_text()); source=ROOT/cfg["selection_source"]
    if digest(source)!=cfg["selection_source_sha256"]: raise RuntimeError("selection source hash mismatch")
    out=ROOT/"outputs"/"v24"/cfg["run_id"]
    if out.exists(): raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True); jobs=[(cfg,s,m,r) for s in cfg["states"] for m in cfg["strategies"] for r in range(cfg["replicates_per_condition"])]
    rows=[]; partial=out/"V24_DESIGN_REPLICATES.partial.csv"
    with open(partial,"w",newline="",encoding="utf-8") as f, ProcessPoolExecutor(max_workers=cfg["worker_processes"]) as pool:
        w=None
        for i,future in enumerate(as_completed([pool.submit(task,j) for j in jobs]),1):
            row=future.result()
            if w is None: w=csv.DictWriter(f,fieldnames=list(row)); w.writeheader()
            w.writerow(row); f.flush(); rows.append(row)
            if i%100==0: print(f"completed {i}/{len(jobs)}",flush=True)
    rows.sort(key=lambda r:(r["state_id"],r["strategy"],r["replicate"])); final=out/"V24_DESIGN_REPLICATES.csv"; write_csv(final,rows); partial.unlink()
    grouped=defaultdict(list)
    for r in rows: grouped[(r["state_id"],r["strategy"])].append(r)
    summaries=[]
    for (sid,strategy),g in sorted(grouped.items()):
        e=np.array([r["error_riu"] for r in g]); mse=float(np.mean(e*e)); state=next(s for s in cfg["states"] if s["state_id"]==sid); fn,p,_,steps,bins=build_model(cfg,state,strategy); _,cols=finite_difference_columns(fn,p,steps); info=information_diagnostics(cols)
        summaries.append({"state_id":sid,"strategy":strategy,"n":len(g),"bias_riu":float(np.mean(e)),"bias_mcse_riu":float(np.std(e,ddof=1)/math.sqrt(len(e))),"mse_riu2":mse,"rmse_riu":math.sqrt(mse),"local_sd_riu":info["local_sd"],"rmse_over_local_sd":math.sqrt(mse)/info["local_sd"],"rank":info["rank"],"parameter_count":info["parameter_count"],"condition_number_whitened_jacobian":info["condition_number_whitened_jacobian"],"failure_rate":float(np.mean([r["optimizer_failure"] for r in g])),"target_boundary_rate":float(np.mean([r["target_boundary"] for r in g])),"multi_solution_rate":float(np.mean([r["multi_solution"] for r in g]))})
    write_csv(out/"V24_DESIGN_CONFIRMATION_SUMMARY.csv",summaries)
    summary={"status":"V24_DESIGN_CONFIRMATION_COMPLETE","run_id":cfg["run_id"],"seed":cfg["seed"],"task_count":len(jobs),"record_count":len(rows),"expected_record_count":4*4*cfg["replicates_per_condition"],"optimizer_failure_count":sum(r["optimizer_failure"] for r in rows),"target_boundary_count":sum(r["target_boundary"] for r in rows),"multi_solution_count":sum(r["multi_solution"] for r in rows),"elapsed_seconds":time.time()-started,"config_sha256":digest(cp),"selection_source_sha256":digest(source),"records_sha256":digest(final),"interpretation_limit":cfg["interpretation_limit"]}
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n"); print(json.dumps(summary,indent=2),flush=True)

if __name__=="__main__":
    if len(sys.argv)!=2: raise SystemExit("usage: run_v24_design_confirmation.py CONFIG")
    main(sys.argv[1])
