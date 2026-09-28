"""Planar multilayer optics with two independent implementations.

Convention: vacuum wavelength, angle from the interface normal, exp(-i omega t),
and passive materials represented by n + i*k with k >= 0.
"""
from __future__ import annotations

import math
from pathlib import Path
import numpy as np


def nk_bk7(wavelength_nm):
    """SCHOTT N-BK7 Sellmeier model; wavelength is in nm."""
    lam_um = np.asarray(wavelength_nm, dtype=float) / 1000.0
    l2 = lam_um * lam_um
    b1, b2, b3 = 1.03961212, 0.231792344, 1.01046945
    c1, c2, c3 = 0.00600069867, 0.0200179144, 103.560653
    n2 = 1.0 + b1*l2/(l2-c1) + b2*l2/(l2-c2) + b3*l2/(l2-c3)
    return np.sqrt(n2).astype(complex)


def load_johnson_gold(path):
    """Parse the tabulated-nk block in the CC0 refractiveindex.info YAML file."""
    rows = []
    in_data = False
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.strip() == "data: |":
            in_data = True
            continue
        if in_data:
            parts = line.split()
            if len(parts) != 3:
                if rows:
                    break
                continue
            try:
                rows.append(tuple(float(x) for x in parts))
            except ValueError:
                if rows:
                    break
    if len(rows) < 3:
        raise ValueError(f"Could not parse Johnson Au data from {path}")
    arr = np.asarray(rows, dtype=float)
    return arr[:, 0] * 1000.0, arr[:, 1], arr[:, 2]


def nk_gold_johnson(wavelength_nm, source_path):
    wl, n, k = load_johnson_gold(source_path)
    x = np.asarray(wavelength_nm, dtype=float)
    if np.any(x < wl.min()) or np.any(x > wl.max()):
        raise ValueError(f"Johnson Au data only support {wl.min():g}-{wl.max():g} nm")
    return np.interp(x, wl, n) + 1j*np.interp(x, wl, k)


def _positive_kz(n, beta):
    kz = np.sqrt(np.asarray(n, dtype=complex)**2 - np.asarray(beta, dtype=complex)**2 + 0j)
    flip = (np.imag(kz) < -1e-14) | ((np.abs(np.imag(kz)) <= 1e-14) & (np.real(kz) < 0))
    return np.where(flip, -kz, kz)


def _admittance(n, beta, pol):
    kz = _positive_kz(n, beta)
    if pol == "s":
        return kz
    if pol == "p":
        return kz / (np.asarray(n, dtype=complex)**2)
    raise ValueError("pol must be 's' or 'p'")


def interface_r(n_i, n_j, beta, pol="p"):
    yi, yj = _admittance(n_i, beta, pol), _admittance(n_j, beta, pol)
    return (yi - yj) / (yi + yj)


def recursive_reflection(n_list, d_nm, wavelength_nm, angle_deg, pol="p"):
    """Stable Fresnel/scattering recursion; finite layer thicknesses are d_nm."""
    n = np.asarray(n_list, dtype=complex)
    if len(n) != len(d_nm) + 2:
        raise ValueError("len(n_list) must equal len(d_nm)+2")
    angle = np.deg2rad(np.asarray(angle_deg, dtype=float))
    beta = n[0] * np.sin(angle)
    r_eff = interface_r(n[-2], n[-1], beta, pol)
    k0 = 2.0 * np.pi / float(wavelength_nm)
    for j in range(len(n)-2, 0, -1):
        phase = np.exp(2j * k0 * _positive_kz(n[j], beta) * float(d_nm[j-1]))
        rij = interface_r(n[j-1], n[j], beta, pol)
        r_eff = (rij + r_eff * phase) / (1.0 + rij * r_eff * phase)
    return r_eff


def characteristic_reflection(n_list, d_nm, wavelength_nm, angle_deg, pol="p"):
    """Independent characteristic-matrix implementation for verification."""
    n = np.asarray(n_list, dtype=complex)
    if len(n) != len(d_nm) + 2:
        raise ValueError("len(n_list) must equal len(d_nm)+2")
    angles = np.deg2rad(np.atleast_1d(np.asarray(angle_deg, dtype=float)))
    out = np.empty(angles.size, dtype=complex)
    k0 = 2.0 * np.pi / float(wavelength_nm)
    for idx, angle in enumerate(angles):
        beta = n[0] * math.sin(float(angle))
        m = np.eye(2, dtype=complex)
        for j, thickness in enumerate(d_nm, start=1):
            y = _admittance(n[j], beta, pol)
            delta = k0 * _positive_kz(n[j], beta) * float(thickness)
            layer = np.array([[np.cos(delta), -1j*np.sin(delta)/y],
                              [-1j*y*np.sin(delta), np.cos(delta)]], dtype=complex)
            m = m @ layer
        y0, ys = _admittance(n[0], beta, pol), _admittance(n[-1], beta, pol)
        b = m[0,0] + m[0,1]*ys
        c = m[1,0] + m[1,1]*ys
        out[idx] = (y0*b - c) / (y0*b + c)
    return out[0] if np.ndim(angle_deg) == 0 else out


def reflectance_kretschmann(angle_deg, wavelength_nm, thickness_nm, analyte_n, gold_source, pol="p", method="recursive"):
    n_prism = complex(nk_bk7(wavelength_nm))
    n_gold = complex(nk_gold_johnson(wavelength_nm, gold_source))
    n_list = [n_prism, n_gold, complex(analyte_n)]
    solver = recursive_reflection if method == "recursive" else characteristic_reflection
    r = solver(n_list, [thickness_nm], wavelength_nm, angle_deg, pol)
    return np.abs(r)**2

