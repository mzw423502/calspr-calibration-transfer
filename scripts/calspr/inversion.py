"""Bounded multistart Poisson inversion without SciPy.

The deterministic Halton start library is defined only by parameter bounds. It
does not use simulated truth, fulfilling the no-truth-initialization contract.
"""
from __future__ import annotations

import numpy as np


def poisson_deviance_half(y, mu):
    mu = np.asarray(mu, dtype=float)
    y = np.asarray(y, dtype=float)
    if np.any(~np.isfinite(mu)) or np.any(mu <= 0):
        return float("inf")
    positive = y > 0
    value = np.sum(mu[~positive])
    value += np.sum(mu[positive] - y[positive] + y[positive]*np.log(y[positive]/mu[positive]))
    return float(value)


def _van_der_corput(index, base):
    value, denom = 0.0, 1.0
    while index:
        index, rem = divmod(index, base)
        denom *= base
        value += rem/denom
    return value


def halton_starts(dim, count):
    primes = [2, 3, 5, 7, 11, 13, 17, 19]
    if dim > len(primes):
        raise ValueError("increase prime library")
    starts = [np.full(dim, 0.5)]
    for i in range(1, count):
        starts.append(np.array([_van_der_corput(i, primes[j]) for j in range(dim)]))
    return starts


def bounded_pattern_fit(objective_physical, bounds, settings):
    bounds = np.asarray(bounds, dtype=float)
    lower, width = bounds[:, 0], bounds[:, 1]-bounds[:, 0]
    if np.any(width <= 0):
        raise ValueError("invalid bounds")

    def objective_unit(u):
        return objective_physical(lower + np.clip(u, 0, 1)*width)

    solutions = []
    for start_id, start in enumerate(halton_starts(len(bounds), settings["halton_starts"])):
        u = np.asarray(start, dtype=float)
        best = objective_unit(u)
        step = float(settings["initial_step_fraction"])
        evaluations = 1
        for _ in range(settings["max_levels"]):
            any_improvement = False
            for _pass in range(settings["passes_per_level"]):
                improved_this_pass = False
                for j in range(len(u)):
                    candidates = []
                    for direction in (-1.0, 1.0):
                        trial = u.copy()
                        trial[j] = np.clip(trial[j] + direction*step, 0.0, 1.0)
                        value = objective_unit(trial)
                        evaluations += 1
                        candidates.append((value, trial))
                    value, trial = min(candidates, key=lambda z: z[0])
                    if value < best - 1e-10:
                        u, best = trial, value
                        improved_this_pass = True
                        any_improvement = True
                if not improved_this_pass:
                    break
            step *= float(settings["step_shrink"])
        solutions.append({"start_id": start_id, "params": lower+u*width,
                          "objective": best, "final_step": step, "evaluations": evaluations,
                          "improved": any_improvement})
    solutions.sort(key=lambda z: z["objective"])
    return solutions


def fisher_scoring_refine(y, model_physical, initial, bounds, max_iter=40, derivative_step=2e-5):
    """Refine one bounded solution with damped Fisher scoring in unit coordinates."""
    bounds = np.asarray(bounds, dtype=float)
    lower, width = bounds[:, 0], bounds[:, 1]-bounds[:, 0]
    u = np.clip((np.asarray(initial, dtype=float)-lower)/width, 0.0, 1.0)

    def model_unit(v):
        return np.asarray(model_physical(lower+np.clip(v, 0, 1)*width), dtype=float)

    mu = model_unit(u)
    best = poisson_deviance_half(y, mu)
    damping = 1e-6
    accepted = 0
    for iteration in range(max_iter):
        jac = np.empty((len(mu), len(u)), dtype=float)
        for j in range(len(u)):
            plus, minus = u.copy(), u.copy()
            plus[j] = min(1.0, plus[j]+derivative_step)
            minus[j] = max(0.0, minus[j]-derivative_step)
            denom = plus[j]-minus[j]
            jac[:, j] = (model_unit(plus)-model_unit(minus))/denom
        score = jac.T @ ((np.asarray(y)-mu)/mu)
        fisher = jac.T @ (jac/mu[:, None])
        scale = np.maximum(np.diag(fisher), 1.0)
        try:
            step = np.linalg.solve(fisher + damping*np.diag(scale), score)
        except np.linalg.LinAlgError:
            damping *= 100.0
            continue
        improved = False
        alpha = 1.0
        for _ in range(16):
            trial = np.clip(u + alpha*step, 0.0, 1.0)
            trial_mu = model_unit(trial)
            value = poisson_deviance_half(y, trial_mu)
            if value < best - 1e-10:
                u, mu, best = trial, trial_mu, value
                accepted += 1
                damping = max(damping/3.0, 1e-12)
                improved = True
                break
            alpha *= 0.5
        if not improved:
            damping *= 10.0
        if np.max(np.abs(alpha*step)) < 2e-8:
            break
    return {"params": lower+u*width, "objective": best, "iterations": iteration+1,
            "accepted_steps": accepted, "damping": damping}
