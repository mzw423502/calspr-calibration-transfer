#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from calspr.optics import (characteristic_reflection, interface_r, nk_bk7,
                           nk_gold_johnson, recursive_reflection,
                           reflectance_kretschmann, _admittance)


def check(name, condition, detail=""):
    if not bool(condition):
        raise AssertionError(f"{name}: {detail}")
    print(f"PASS {name} {detail}")


def main():
    gold = ROOT / "data/raw/materials/Au_Johnson.yml"
    check("BK7_633_range", 1.50 < nk_bk7(633).real < 1.53, str(nk_bk7(633)))
    au = nk_gold_johnson(633, gold)
    check("Au_passive", au.imag > 0, str(au))

    for pol in ("s", "p"):
        beta = 1.0*math.sin(math.radians(25))
        r_direct = interface_r(1.0, 1.5, beta, pol)
        r_rec = recursive_reflection([1.0, 1.5], [], 633, 25, pol)
        check(f"single_interface_{pol}", abs(r_direct-r_rec) < 1e-14, str(abs(r_direct-r_rec)))

    rs = recursive_reflection([1.0, 1.5], [], 633, 0, "s")
    rp = recursive_reflection([1.0, 1.5], [], 633, 0, "p")
    check("normal_incidence_ps_reflectance", abs(abs(rs)**2-abs(rp)**2) < 1e-14)

    beta = 1.0*math.sin(math.radians(30))
    for pol in ("s", "p"):
        yi, yj = _admittance(1.0, beta, pol), _admittance(1.5, beta, pol)
        r = (yi-yj)/(yi+yj)
        t = 2*yi/(yi+yj)
        balance = abs(r)**2 + (yj/yi).real*abs(t)**2
        check(f"lossless_flux_{pol}", abs(balance-1) < 1e-12, str(balance))

    angles = np.linspace(60, 80, 101)
    n_stack = [complex(nk_bk7(633)), complex(au), 1.333+0j]
    rr = recursive_reflection(n_stack, [50.0], 633, angles, "p")
    rc = characteristic_reflection(n_stack, [50.0], 633, angles, "p")
    check("independent_solver_reflectance", np.max(np.abs(np.abs(rr)**2-np.abs(rc)**2)) < 1e-10,
          str(np.max(np.abs(np.abs(rr)**2-np.abs(rc)**2))))

    r0 = recursive_reflection(n_stack, [0.0], 633, angles, "p")
    rd = recursive_reflection([n_stack[0], n_stack[2]], [], 633, angles, "p")
    check("zero_thickness_degeneracy", np.max(np.abs(r0-rd)) < 1e-10, str(np.max(np.abs(r0-rd))))

    rthick = recursive_reflection(n_stack, [5000.0], 633, angles, "p")
    rmetal = recursive_reflection(n_stack[:2], [], 633, angles, "p")
    check("thick_metal_limit", np.max(np.abs(np.abs(rthick)**2-np.abs(rmetal)**2)) < 1e-9)

    refl = reflectance_kretschmann(angles, 633, 50, 1.333, gold)
    check("reflectance_physical", np.all((refl >= 0) & (refl <= 1+1e-10)), f"min={refl.min()} max={refl.max()}")
    check("spr_feature_present", refl.min() < 0.5, f"minimum={refl.min()}")

    result = {"tests": 11, "status": "PASS", "numpy": np.__version__}
    (ROOT / "outputs" / "test_summary.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()

