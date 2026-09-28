#!/usr/bin/env python3
from __future__ import annotations
import json, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if len(sys.argv)!=3:raise SystemExit("usage: prepare_v24_design_smoke.py SOURCE TARGET")
source=ROOT/sys.argv[1]; target=ROOT/sys.argv[2]
if target.exists():raise FileExistsError(f"refusing to overwrite {target}")
cfg=json.loads(source.read_text()); cfg["run_id"]="v24_design_confirmation_smoke_20260914_s1"; cfg["run_class"]="v24_design_confirmation_smoke"; cfg["replicates_per_condition"]=1; cfg["worker_processes"]=1; cfg["states"]=[cfg["states"][0]]; cfg["interpretation_limit"]="Smoke and runtime check only; not inferential evidence."
target.write_text(json.dumps(cfg,indent=2)+"\n"); print(target)
