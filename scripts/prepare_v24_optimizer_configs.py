#!/usr/bin/env python3
"""Prepare paired robust/fast optimizer validation and final V2.4 configs."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_new(path, obj):
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    path.write_text(json.dumps(obj, indent=2)+"\n", encoding="utf-8")
    return {"path": str(path), "sha256": digest(path)}


def main(source_config):
    source = Path(source_config)
    if not source.is_absolute():
        source = ROOT / source
    base = json.loads(source.read_text(encoding="utf-8"))
    outputs = []

    robust = json.loads(json.dumps(base))
    robust["run_id"] = "v24_optimizer_validation_robust_20260914_v1"
    robust["run_class"] = "optimizer_validation_robust_reference"
    robust["replicates_per_condition"] = 1
    robust["worker_processes"] = 4
    robust["interpretation_limit"] = "Optimizer validation only; not inferential evidence."
    outputs.append(write_new(ROOT/"configs"/"v24_optimizer_validation_robust.json", robust))

    fast = json.loads(json.dumps(robust))
    fast["run_id"] = "v24_optimizer_validation_fast_20260914_v1"
    fast["run_class"] = "optimizer_validation_candidate_fast"
    fast["fit"]["halton_starts"] = 1
    fast["fit"]["passes_per_level"] = 1
    fast["fit"]["max_levels"] = 9
    outputs.append(write_new(ROOT/"configs"/"v24_optimizer_validation_fast.json", fast))

    final = json.loads(json.dumps(base))
    final["run_id"] = "v24_confirmation_20260914_c1_fast_validated"
    final["run_class"] = "frozen_v24_paired_confirmation_fast_after_validation"
    final["fit"]["halton_starts"] = 1
    final["fit"]["passes_per_level"] = 1
    final["fit"]["max_levels"] = 9
    final["optimizer_validation_gate"] = {
        "robust_config": "configs/v24_optimizer_validation_robust.json",
        "fast_config": "configs/v24_optimizer_validation_fast.json",
        "required_before_launch": True,
        "acceptance": "all paired q estimates differ by <=1e-7 RIU; no objective is worse by >1e-6 half-deviance; zero failures, boundary hits, and multi-solution flags"
    }
    outputs.append(write_new(ROOT/"configs"/"v24_confirmation_frozen_fast_v1.json", final))
    print(json.dumps({"source": str(source), "source_sha256": digest(source),
                      "outputs": outputs}, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: prepare_v24_optimizer_configs.py SOURCE_FROZEN_CONFIG")
    main(sys.argv[1])
