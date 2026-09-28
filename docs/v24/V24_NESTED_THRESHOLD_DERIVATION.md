# CalSPR V2.4 nested-model threshold identity

## Scope

This note verifies a standard block-matrix identity for the actual CalSPR T0/T1 observation model. It is not a new statistical theorem. It applies only when T0 and T1 use the same observations and constraints, the additional mismatch parameter is scalar, the model is locally correct and differentiable, the Fisher matrix is full rank, and the operating point is regular and interior.

## Actual Fisher partition

Let `theta = (q, eta)` contain the target refractive index and every nuisance parameter fitted by rigid transfer T0. Let `delta` be the one target-minus-source thickness difference additionally fitted by T1. For the same stacked source-reference and target-sample Poisson observations, partition the full T1 Fisher matrix at `delta = 0` as

`F_full = [[A, b], [b^T, c]]`.

Here `A` is exactly the T0 Fisher matrix, `b` contains cross-information with the additional thickness difference, and `c` is its unprofiled information. The profiled information for `delta` is the Schur complement

`S_delta = c - b^T A^{-1} b`.

With `e_q` selecting the target coordinate, the actual sign convention used by the generating and fitting code gives

`K_delta = e_q^T A^{-1} b`.

For a small true positive target-minus-source thickness difference omitted by T0, this coefficient gives the first-order T0 target shift `K_delta delta`. The legacy point-response calculation reproduced the positive V2.2 coefficient, confirming the sign convention.

## Variance and threshold identities

The T0 local target variance is

`v0 = (A^{-1})_qq`.

Block inversion of `F_full` gives the T1 target variance `v1` and

`v1 - v0 = K_delta^2 / S_delta`.

Under the local risk approximations

`MSE_T0(delta) = v0 + (K_delta delta)^2`,

`MSE_T1(delta) = v1`,

the nonzero crossing is

`delta_star = sqrt(v1-v0)/abs(K_delta)`.

Substitution of the block-inverse identity yields

`delta_star = 1/sqrt(S_delta)`.

The second form shows that, under the stated regularity conditions, the local mismatch scale is the Cramer-Rao standard deviation for estimating the added thickness-difference coordinate in T1. It remains design-specific and must not be called a measured fabrication tolerance.

## Numerical verification

The main V2.4 response is a 0.10-degree top-hat angular average over the 201 historical centres. Nine anchors cross nominal thicknesses 45, 50, and 55 nm with total incident-photon budgets 1e5, 1e6, and 1e7. All five-column whitened Jacobians were full rank. Across all anchors, both the variance-increment and threshold identities agreed to much better than `1e-6` relative tolerance; observed errors were of order `1e-15`.

Predicted thresholds (nm) were:

| Thickness (nm) | 1e5 photons | 1e6 photons | 1e7 photons |
|---:|---:|---:|---:|
| 45 | 0.42910 | 0.13569 | 0.04291 |
| 50 | 0.31025 | 0.09811 | 0.03102 |
| 55 | 0.26722 | 0.08450 | 0.02672 |

Within each thickness, `S_delta` scales proportionally with photon count and the predicted threshold follows the expected inverse-square-root photon scaling. The thickness dependence is not monotone in sensitivity alone; it reflects the complete target/nuisance geometry of the stated observation model.

## Pre-noise confirmation design

For every anchor, noisy paired inversion uses

`delta/delta_pred = 0, +/-0.5, +/-1, +/-2`.

T0 and T1 are fitted to the same Poisson realization in each replicate. The primary paired quantity is

`D = squared_error_T0 - squared_error_T1`,

so positive D favors T1. The local prediction is

`D_pred = (v1-v0)[(delta/delta_pred)^2 - 1]`.

The initial confirmation contains 100 replicates per condition, giving 6,300 paired tasks and 12,600 fits. An extension to 200 replicates is permitted only for the whole nine-anchor group if a predeclared precision rule is triggered at absolute normalized mismatch 0.5 or 2. No extension is allowed merely to force resolution at normalized mismatch 1, where the local predicted risk difference is zero.

## Failure conditions

No threshold is interpreted when the full Fisher matrix is rank deficient, `S_delta <= 0`, the fit is controlled by a parameter bound, multiple solutions remain unresolved, the generating and fitted observation models differ, or T0 and T1 do not use the same observations and constraints. Bare-prism P has no gold-thickness difference and therefore no CalSPR thickness-transfer threshold.

