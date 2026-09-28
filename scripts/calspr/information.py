"""Poisson observation and nuisance-projected information calculations."""
from __future__ import annotations

import numpy as np
from .optics import reflectance_kretschmann


def means(angles, wavelength_nm, d_nm, q, gold_source, photons, background, angle_zero=0.0, log_gain=0.0):
    r = reflectance_kretschmann(np.asarray(angles)+angle_zero, wavelength_nm, d_nm, q, gold_source)
    per_bin = float(photons) / len(angles)
    return per_bin * (np.exp(log_gain)*r + float(background))


def centered_derivative(fn, x, h):
    return (fn(x+h)-fn(x-h))/(2.0*h)


def effective_information(derivatives, target=0, prior_rows=None, rtol=1e-10):
    """Project the whitened target derivative away from the nuisance column space."""
    a = np.asarray(derivatives, dtype=float)
    if prior_rows is not None:
        a = np.vstack([a, np.asarray(prior_rows, dtype=float)])
    g = a[:, target]
    nuisance = np.delete(a, target, axis=1)
    if nuisance.shape[1] == 0:
        return float(g@g), 0, 1.0
    u, s, _ = np.linalg.svd(nuisance, full_matrices=False)
    tol = rtol * (s[0] if len(s) and s[0] else 1.0)
    rank = int(np.sum(s > tol))
    resid = g - u[:, :rank] @ (u[:, :rank].T @ g)
    cos_resid = float(np.linalg.norm(resid) / np.linalg.norm(g)) if np.linalg.norm(g) else 0.0
    return float(resid@resid), rank, cos_resid


def scenario_derivatives(cfg, d_nm, regime, gold_source):
    angles = np.linspace(cfg["angle_min_deg"], cfg["angle_max_deg"], cfg["angle_bins"])
    w = cfg["wavelength_nm"]
    q0 = cfg["baseline_index"]
    dq = cfg["sample_delta_index"]
    total = cfg["total_incident_photons"]
    bg = cfg["background_fraction"]
    h = cfg["finite_difference_steps"]

    if regime == "absolute":
        def vector(p):
            q, d, a0, lg = p
            return means(angles, w, d, q, gold_source, total, bg, a0, lg)
        p = np.array([q0+dq, d_nm, 0.0, 0.0])
        steps = np.array([h["refractive_index"], h["thickness_nm"], h["angle_zero_deg"], h["log_gain"]])
        labels = ["q", "d_nm", "angle_zero_deg", "log_gain"]
    elif regime == "same_device_differential":
        def vector(p):
            delta, d, a0, lg = p
            base = means(angles, w, d, q0, gold_source, total/2, bg, a0, lg)
            sample = means(angles, w, d, q0+delta, gold_source, total/2, bg, a0, lg)
            return np.concatenate([base, sample])
        p = np.array([dq, d_nm, 0.0, 0.0])
        steps = np.array([h["refractive_index"], h["thickness_nm"], h["angle_zero_deg"], h["log_gain"]])
        labels = ["delta_q", "shared_d_nm", "shared_angle_zero_deg", "shared_log_gain"]
    elif regime == "cross_device_transfer":
        def vector(p):
            q, d_cal, mismatch, a0, lg = p
            cal = means(angles, w, d_cal, q0, gold_source, total/2, bg, a0, lg)
            sample = means(angles, w, d_cal+mismatch, q, gold_source, total/2, bg, a0, lg)
            return np.concatenate([cal, sample])
        p = np.array([q0+dq, d_nm, 0.0, 0.0, 0.0])
        steps = np.array([h["refractive_index"], h["thickness_nm"], h["thickness_nm"], h["angle_zero_deg"], h["log_gain"]])
        labels = ["q", "calibration_d_nm", "device_mismatch_nm", "shared_angle_zero_deg", "shared_log_gain"]
    else:
        raise ValueError(regime)

    mu = vector(p)
    cols = []
    for j, step in enumerate(steps):
        cols.append(centered_derivative(lambda x: vector(np.r_[p[:j], x, p[j+1:]]), p[j], step) / np.sqrt(mu))
    return np.column_stack(cols), labels, angles, mu


def calibration_prior_rows(cfg, labels):
    sig = cfg["calibration_scenarios"]
    rows = []
    mapping = {
        "d_nm": sig["thickness_sigma_nm"],
        "shared_d_nm": sig["thickness_sigma_nm"],
        "calibration_d_nm": sig["thickness_sigma_nm"],
        "device_mismatch_nm": sig["transfer_mismatch_sigma_nm"],
        "angle_zero_deg": sig["angle_zero_sigma_deg"],
        "shared_angle_zero_deg": sig["angle_zero_sigma_deg"],
        "log_gain": sig["log_gain_sigma"],
        "shared_log_gain": sig["log_gain_sigma"]
    }
    for name, sigma in mapping.items():
        if name in labels:
            row = np.zeros(len(labels))
            row[labels.index(name)] = 1.0/float(sigma)
            rows.append(row)
    return np.asarray(rows)

