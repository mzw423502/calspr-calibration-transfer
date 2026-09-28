#!/usr/bin/env python3
"""Create one immutable retry config after an interrupted launcher attempt."""
from __future__ import annotations
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if len(sys.argv) != 4:
    raise SystemExit("usage: prepare_v24_retry_config.py SOURCE TARGET RUN_ID")
source = ROOT / sys.argv[1]
target = ROOT / sys.argv[2]
if target.exists():
    raise FileExistsError(f"refusing to overwrite {target}")
cfg = json.loads(source.read_text(encoding="utf-8"))
cfg["run_id"] = sys.argv[3]
cfg["retry_of_interrupted_run_id"] = "v24_optimizer_validation_robust_20260914_v1"
cfg["retry_reason"] = "The first launcher shell was interrupted while correcting a PID-path error; partial outputs are preserved and excluded."
target.write_text(json.dumps(cfg, indent=2)+"\n", encoding="utf-8")
print(target)
