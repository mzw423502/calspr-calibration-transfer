#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()


def main(config_path):
    cp=Path(config_path).resolve(); cfg=json.loads(cp.read_text(encoding="utf-8"))
    src=ROOT/cfg["input"]; out=ROOT/cfg["output_directory"]
    if out.exists(): raise FileExistsError(f"refusing to overwrite {out}")
    out.mkdir(parents=True)
    with open(src,encoding="utf-8") as f: rows=list(csv.DictReader(f))
    corrected=[]
    for r in rows:
        z=dict(r)
        for prefix in ("observation","scenario_prior"):
            old=f"{prefix}_separation_fraction"
            value=float(z.pop(old))
            z[f"{prefix}_target_derivative_norm_retained"]=value
            z[f"{prefix}_target_information_fraction_retained"]=value*value
        corrected.append(z)
    dest=out/"information_scan_corrected.csv"
    with open(dest,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(corrected[0])); w.writeheader(); w.writerows(corrected)
    report={
        "status":"COMPLETE","run_id":cfg["run_id"],"source":str(src),"source_sha256":sha(src),
        "output":str(dest),"output_sha256":sha(dest),"rows":len(corrected),
        "transformation":"Renamed stored separation_fraction values as derivative-norm ratios and added their squares as retained information fractions. No historical file or other numeric column was changed."
    }
    (out/"summary.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    if len(sys.argv)!=2: raise SystemExit("usage: run_correct_information_labels.py CONFIG")
    main(sys.argv[1])
