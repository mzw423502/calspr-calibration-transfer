# Local bias-variance derivation for calibration transfer

## Setup

Let the stacked raw Poisson counts have nominal mean `mu(q, eta, delta)`. The target is refractive index `q`; `eta` contains nuisance parameters still fitted by the rigid transfer model; `delta` is a target-minus-source device difference omitted when the rigid model fixes `delta = 0`.

At a regular interior point with `mu0 > 0`, define

`W = diag(1 / mu0)`,

`u = sqrt(W) dmu/dq`, `V = sqrt(W) dmu/deta`, and `v = sqrt(W) dmu/ddelta`.

Let `P_V` be the orthogonal projector onto the column space of `V`, and define the nuisance-residualized derivatives

`r_q = (I - P_V) u`,

`r_delta = (I - P_V) v`.

The effective Fisher information for `q` in the rigid model is

`I_rigid = r_q^T r_q`.

The implementation obtains the projector from an SVD after scaling nuisance columns to unit norm. Column scaling does not change their span but improves numerical conditioning. Rank and finite-difference convergence must still be checked.

## First-order omitted-mismatch bias

For a small true mismatch `delta`, linearize the true mean around zero. After profiling the fitted nuisance directions, the residual mean perturbation is approximately `r_delta delta`. The local score equation for a shift `Delta q` is therefore

`r_q^T (r_delta delta - r_q Delta q) approximately 0`.

Solving gives

`Bias(q_hat) approximately K_delta delta`,

where

`K_delta = (r_q^T r_delta) / (r_q^T r_q)`.

This is standard local omitted-variable/estimating-equation theory applied to the stated Poisson model, not a new statistical theorem. CalSPR tests its usefulness as a pre-simulation predictor. `K_delta` is computed from forward-model derivatives at zero mismatch, then evaluated against deterministic pseudo-true fits at independent positive and negative mismatch values. It is not estimated by regressing the already observed bias curve.

For finite mismatch, the pseudo-true parameter minimizes the Poisson Kullback-Leibler discrepancy, equivalently the Poisson deviance with the true mean used as deterministic data. Those fits describe nonlinear model-misspecification bias and are not stochastic replicates.

## Variance cost of estimating the difference

When `delta` is explicitly fitted without new observations or an informative prior, project `r_q` again against the effective mismatch direction. With one additional effective nuisance direction, let

`c = (r_q^T r_delta) / (||r_q|| ||r_delta||)`.

Then

`I_joint = ||r_q||^2 (1 - c^2)`.

The corresponding local variances are `v_rigid = 1/I_rigid` and `v_joint = 1/I_joint`. Under model match, interior regularity, the same observation design, and locally negligible estimator bias,

`MSE_rigid approximately v_rigid + (K_delta delta)^2`,

`MSE_joint approximately v_joint`.

If `v_joint >= v_rigid` and `K_delta` is nonzero, the approximate strategy-crossing scale is

`abs(delta_star) approximately sqrt(v_joint - v_rigid) / abs(K_delta)`.

No threshold is reported when the numerator is negative, `K_delta` is effectively zero, the information is singular, or the fit is controlled by bounds or multiple solutions. A target-device reference (T2) changes the observation design while retaining the same total photon budget, so its risk comparison is computed separately; it is not inserted into the one-added-direction identity without re-derivation.

## Added target-angle difference

For chip reinstallation, a target-specific angle offset `tau` is treated analogously. A rigid model can have first-order coefficient

`K_tau = (r_q^T r_tau) / (r_q^T r_q)`.

T1-angle and T2-angle explicitly fit both thickness difference and target angle offset. Their loss of information is reported directly from the full projected derivative matrix. This scenario is distinct from the stable-common-angle scenario; neither is used as evidence that the other is false.

## Applicability limits

The formulas are local and conditional on the selected material model, angle window, binning, photon allocation, background definition, fitted parameter set, and regular interior point. They do not by themselves cover strong nonlinearity, boundary estimators, multiple optima, material-model mismatch, finite-reference uncertainty, or empirical device populations. Monte Carlo confirmation must use conditions and seeds not used to select the strategy.

