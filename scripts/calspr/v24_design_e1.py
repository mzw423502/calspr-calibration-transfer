"""CalSPR V2.4 angular-response and information-design utilities.

The legacy model evaluates reflectance at angle-bin centres.  This module keeps
that model exactly available as ``point`` and adds explicitly labelled finite
angular-response scenarios based on normalized top-hat integration.  Nothing
in this file assigns an empirical instrument response.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .optics import (
    interface_r,
    nk_bk7,
    nk_gold_johnson,
    reflectance_kretschmann,
)


@dataclass(frozen=True)
class ResponseScenario:
    name: str
    kind: str
    width_deg: float
    quadrature_order: int = 12


def angle_centres(cfg, phase_fraction=0.0):
    """Return the historical 201 centres, optionally shifted for phase checks."""
    centres = np.linspace(cfg["angle_min_deg"], cfg["angle_max_deg"], cfg["angle_bins"])
    spacing = float(centres[1] - centres[0])
    return centres + float(phase_fraction) * spacing


def _top_hat_average(fn, centres, width_deg, order):
    if width_deg <= 0:
        return np.asarray(fn(centres), dtype=float)
    nodes, weights = np.polynomial.legendre.leggauss(int(order))
    sample_angles = centres[:, None] + 0.5 * float(width_deg) * nodes[None, :]
    values = np.asarray(fn(sample_angles.reshape(-1)), dtype=float).reshape(len(centres), -1)
    return 0.5 * (values @ weights)


def gold_reflectance(centres, wavelength_nm, d_nm, q, gold_source, scenario):
    fn = lambda a: reflectance_kretschmann(a, wavelength_nm, d_nm, q, gold_source)
    if scenario.kind == "point":
        return np.asarray(fn(centres), dtype=float)
    if scenario.kind == "top_hat":
        return _top_hat_average(fn, np.asarray(centres), scenario.width_deg,
                                scenario.quadrature_order)
    raise ValueError(f"unknown angular response kind: {scenario.kind}")


def bare_prism_reflectance(centres, wavelength_nm, q, scenario):
    n_prism = complex(nk_bk7(wavelength_nm))

    def fn(a):
        beta = n_prism * np.sin(np.deg2rad(np.asarray(a, dtype=float)))
        return np.abs(interface_r(n_prism, complex(q), beta, pol="p")) ** 2

    if scenario.kind == "point":
        return np.asarray(fn(centres), dtype=float)
    if scenario.kind == "top_hat":
        return _top_hat_average(fn, np.asarray(centres), scenario.width_deg,
                                scenario.quadrature_order)
    raise ValueError(f"unknown angular response kind: {scenario.kind}")


def poisson_mean(reflectance, incident_photons, background_fraction, log_gain=0.0):
    reflectance = np.asarray(reflectance, dtype=float)
    if len(reflectance) == 0:
        return np.empty(0, dtype=float)
    per_bin = float(incident_photons) / len(reflectance)
    return per_bin * (np.exp(log_gain) * reflectance + float(background_fraction))


def physical_windows(cfg, centres, gold_source):
    """Freeze disjoint C/S/O masks from the nominal design, never simulated truth.

    C is centred on the nominal prism/reference critical angle.  S is centred
    on the nominal point-response SPR minimum at the planned sample index.  O
    is the complement.  Widths are configuration values and masks cover F.
    """
    centres = np.asarray(centres, dtype=float)
    q0 = float(cfg["baseline_index"])
    q_plan = q0 + float(cfg["sample_delta_index"])
    n_prism = float(np.real(nk_bk7(cfg["wavelength_nm"])))
    critical = math.degrees(math.asin(q0 / n_prism))
    nominal = reflectance_kretschmann(
        centres, cfg["wavelength_nm"], cfg["window_definition"]["nominal_thickness_nm"],
        q_plan, gold_source,
    )
    spr = float(centres[int(np.argmin(nominal))])
    c_half = float(cfg["window_definition"]["critical_half_width_deg"])
    s_half = float(cfg["window_definition"]["spr_half_width_deg"])
    c_mask = np.abs(centres - critical) <= c_half + 1e-12
    s_mask = (np.abs(centres - spr) <= s_half + 1e-12) & ~c_mask
    o_mask = ~(c_mask | s_mask)
    if not np.all(c_mask | s_mask | o_mask):
        raise AssertionError("window masks do not cover the full grid")
    if np.any(c_mask & s_mask) or np.any(c_mask & o_mask) or np.any(s_mask & o_mask):
        raise AssertionError("window masks overlap")
    return {
        "critical_angle_deg": critical,
        "spr_minimum_deg": spr,
        "C": c_mask,
        "S": s_mask,
        "O": o_mask,
        "CS": c_mask | s_mask,
        "F": np.ones(len(centres), dtype=bool),
        "no_C": ~c_mask,
        "no_S": ~s_mask,
    }


def finite_difference_columns(fn, p, steps):
    p = np.asarray(p, dtype=float)
    mu = np.asarray(fn(p), dtype=float)
    if np.any(~np.isfinite(mu)) or np.any(mu <= 0):
        raise ValueError("non-positive or non-finite Poisson mean")
    cols = []
    for j, h in enumerate(steps):
        plus, minus = p.copy(), p.copy()
        plus[j] += float(h)
        minus[j] -= float(h)
        cols.append((np.asarray(fn(plus)) - np.asarray(fn(minus))) /
                    (2.0 * float(h) * np.sqrt(mu)))
    return mu, np.column_stack(cols)


def information_diagnostics(columns, target=0, rtol=1e-10):
    x = np.asarray(columns, dtype=float)
    fisher = x.T @ x
    column_norms = np.linalg.norm(x, axis=0)
    scaled = x / column_norms if np.all(column_norms > 0) else x.copy()
    s = np.linalg.svd(scaled, compute_uv=False)
    tol = rtol * (s[0] if len(s) and s[0] else 1.0)
    rank = int(np.sum(s > tol))
    condition = float(s[0] / s[-1]) if len(s) and s[-1] > 0 else math.inf
    nuisance = np.delete(x, target, axis=1)
    if nuisance.shape[1]:
        u, ns, _ = np.linalg.svd(nuisance, full_matrices=False)
        nrank = int(np.sum(ns > rtol * (ns[0] if len(ns) and ns[0] else 1.0)))
        residual = x[:, target] - u[:, :nrank] @ (u[:, :nrank].T @ x[:, target])
    else:
        nrank = 0
        residual = x[:, target].copy()
    info = float(residual @ residual)
    variance = 1.0 / info if info > 0 and rank == x.shape[1] else math.inf
    return {
        "fisher": fisher,
        "rank": rank,
        "parameter_count": int(x.shape[1]),
        "nuisance_rank": nrank,
        "condition_number_whitened_jacobian": condition,
        "condition_number_physical_fisher": float(np.linalg.cond(fisher)),
        "effective_information": info,
        "local_variance": variance,
        "local_sd": math.sqrt(variance) if math.isfinite(variance) else math.inf,
        "projected_target_score": residual,
    }


def response_scenarios(cfg):
    out = []
    for item in cfg["angular_response_scenarios"]:
        out.append(ResponseScenario(
            name=item["name"], kind=item["kind"], width_deg=float(item.get("width_deg", 0.0)),
            quadrature_order=int(item.get("quadrature_order", 12)),
        ))
    return out
