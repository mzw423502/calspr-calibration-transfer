#!/usr/bin/env python3
"""Audit candidate fast V2.4 fits against the four-start reference fits."""
from __future__ import annotations
import csv, hashlib, json, sys, time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def digest(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()
def load(p):
    with open(p,newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def key(r):return (r["anchor_id"],r["rho_delta_over_pred"],r["replicate"],r["strategy"])
def main(cp):
    started=time.time(); cp=Path(cp).resolve(); cfg=json.loads(cp.read_text()); rp=ROOT/cfg["robust_records"]; fp=ROOT/cfg["fast_records"]
    robust={key(r):r for r in load(rp)}; fast={key(r):r for r in load(fp)}
    if set(robust)!=set(fast): raise RuntimeError("record keys differ")
    rows=[]
    for k in sorted(robust):
        r,f=robust[k],fast[k]; qdiff=abs(float(f["estimated_q"])-float(r["estimated_q"])); excess=float(f["objective_half_deviance"])-float(r["objective_half_deviance"])
        diag=sum(int(f[x]) for x in ("optimizer_failure","target_boundary","multi_solution"))
        passed=qdiff<=cfg["max_abs_q_difference_riu"] and excess<=cfg["max_fast_objective_excess_half_deviance"] and (not cfg["require_zero_diagnostics"] or diag==0)
        rows.append({"anchor_id":k[0],"rho_delta_over_pred":k[1],"replicate":k[2],"strategy":k[3],"robust_estimated_q":r["estimated_q"],"fast_estimated_q":f["estimated_q"],"abs_q_difference_riu":qdiff,"robust_objective_half_deviance":r["objective_half_deviance"],"fast_objective_half_deviance":f["objective_half_deviance"],"fast_objective_excess_half_deviance":excess,"fast_diagnostic_count":diag,"pass":int(passed)})
    out=ROOT/"outputs"/"v24"/cfg["run_id"]
    if out.exists(): raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True); csvp=out/"V24_OPTIMIZER_VALIDATION.csv"
    with open(csvp,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary={"status":"PASS" if all(r["pass"] for r in rows) else "FAIL","run_id":cfg["run_id"],"record_count":len(rows),"failed_record_count":sum(1-r["pass"] for r in rows),"max_abs_q_difference_riu":max(r["abs_q_difference_riu"] for r in rows),"max_fast_objective_excess_half_deviance":max(r["fast_objective_excess_half_deviance"] for r in rows),"fast_diagnostic_count":sum(r["fast_diagnostic_count"] for r in rows),"elapsed_seconds":time.time()-started,"config_sha256":digest(cp),"robust_records_sha256":digest(rp),"fast_records_sha256":digest(fp),"audit_csv_sha256":digest(csvp)}
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n");print(json.dumps(summary,indent=2))
    if summary["status"]!="PASS":raise SystemExit(2)
if __name__=="__main__":
    if len(sys.argv)!=2:raise SystemExit("usage: run_v24_optimizer_audit.py CONFIG")
    main(sys.argv[1])
