# CalSPR V2.4 measurement response and frozen angular windows

## Status and scope

This document records the observation semantics and the pre-noise angular-window definitions used by CalSPR V2.4. The calculations are conditional simulations of an angular-count model. No angular response, drift distribution, fabrication tolerance, or achieved instrumental resolution was measured.

Authoritative deterministic run: `outputs/v24/v24_scan_20260914_e4_final`.

## Historical 201-bin observation

The V2.2 inverse model used 201 equally spaced angle centres from 60 to 80 degrees, inclusive. The spacing is 0.1 degrees. For an exposure assigned `N_e` incident photons, every centre received `N_e/201` incident photons and its mean reflected count was

`mu_j = (N_e/201) [exp(g) R(theta_j + a_0) + b]`,

with `b = 0.001`. Thus the historical model was centre-point sampling, not numerical integration over a physical pixel and not convolution with a measured point-spread function. The exact centre-point implementation is retained as `legacy_point`.

The V2.4 re-run reproduced the V2.2 50-nm, 10-million-photon threshold: `0.0309855306 nm` versus the recorded `0.0309855387 nm` (relative difference `2.62e-7`). This is a regression check, not a new result.

## Named finite-angular-response scenarios

V2.4 adds normalized top-hat angular averages with widths 0.05, 0.10, and 0.20 degrees. The 0.10-degree response is the main scenario because it equals the historical centre spacing; 0.05 and 0.20 degrees are bounded sensitivity endpoints. These are numerical response scenarios, not estimates of a particular instrument.

Each top-hat mean is evaluated by Gauss-Legendre quadrature. Any bin containing the ideal prism/analyte critical angle is split at that square-root branch point before quadrature. The main calculation uses order 16. At the 50-nm, 10-million-photon anchor, predicted thresholds at orders 8, 16, and 32 were `0.0310245839`, `0.0310245852`, and `0.0310245854 nm`, respectively. The corresponding C- and bare-prism information values changed by less than `5e-5` relatively from order 16 to 32.

Centered finite-difference multipliers 0.5, 1, and 2 and grid phases -0.5, 0, and +0.5 of one 0.1-degree spacing were checked. The full-curve threshold and the SPR-window result were stable to grid phase. Critical-angle-only and bare-prism information changed materially with phase, as expected for a sharp ideal critical feature sampled by a coarse grid. That phase dependence is retained as a limitation rather than averaged away.

## Frozen disjoint windows

Windows were defined before noisy outcomes from nominal physical quantities:

- Nominal prism/reference critical angle: `61.6206773 degrees`, calculated from BK7 and `q0 = 1.3330`.
- Nominal SPR minimum: `73.0 degrees`, obtained from the point-response 50-nm Johnson-Christy model at the planned sample index `q = 1.3380`.
- C: centres within 1.0 degree of the nominal critical angle, 60.7-62.6 degrees, 20 bins.
- S: centres within 2.0 degrees of the nominal SPR minimum, excluding C, 71.0-75.0 degrees, 41 bins.
- O: all remaining centres, 140 bins.
- F: all 201 bins.

C, S, and O are mutually exclusive and cover F exactly. No simulated true thickness difference or noisy target estimate enters these definitions.

## Attribution without photon reallocation

For F, no-C, no-S, C, S, CS, and O, photons per retained bin were kept equal to the F design. Photons belonging to removed bins were discarded, not reassigned. All four fitted columns (`q`, thickness, angle zero, and log gain) were projected jointly. Separately profiled regional information values were not added.

At 10 million original incident photons under the main response:

| Design | Fraction of original photons used | Effective information | Local SD (RIU) |
|---|---:|---:|---:|
| F | 1.000 | 3.5805e9 | 1.6712e-5 |
| no C | 0.9005 | 3.5709e9 | 1.6734e-5 |
| no S | 0.7960 | 6.8004e8 | 3.8347e-5 |
| C | 0.0995 | 4.4718e6 | 4.7289e-4 |
| S | 0.2040 | 2.6647e9 | 1.9372e-5 |
| CS | 0.3035 | 2.7123e9 | 1.9201e-5 |
| O | 0.6965 | 6.7325e8 | 3.8540e-5 |

Within the single joint F projection, regional squared-score contributions were 9.03e6 (C), 2.7603e9 (S), and 8.1109e8 (O), summing to the F effective information. These are contributions under the F nuisance projection, not stand-alone regional Fisher informations.

## Equal-total-budget comparison

F, C, S, CS, and bare prism P each received the same 10-million-photon total budget. Reference fractions 0.20, 0.35, 0.50, 0.65, and 0.80 were screened before noise. CS also screened C weights 0, 0.25, 0.50, 0.75, and 1; both endpoints were included. All strategies used the same 0.10-degree top-hat response and the same reference/sample accounting.

Gold designs fitted `q`, gold thickness, angle zero, and log gain. P omitted gold and therefore had no gold-thickness parameter, but still fitted `q`, angle zero, and log gain. A thickness-transfer threshold is not applicable to P.

| Strategy | Selected reference fraction | Selected C weight | Effective information | Local SD (RIU) |
|---|---:|---:|---:|---:|
| F | 0.50 | N/A | 3.5805e9 | 1.6712e-5 |
| C | 0.50 | N/A | 4.4941e7 | 1.4917e-4 |
| S | 0.50 | N/A | 1.3063e10 | 8.7493e-6 |
| CS | 0.50 | 0.00 | 1.3063e10 | 8.7493e-6 |
| P | 0.50 | N/A | 1.3390e9 | 2.7328e-5 |

The CS optimum lies exactly at the S-only endpoint. In this local conditional design, the resonance window carries most target information, whereas the critical window adds little once thickness, angle zero, and gain are fitted. P retains bulk-index information without a gold-film nuisance but does not outperform the selected S window. These statements are local design results, not empirical claims about all SPR or critical-angle instruments.

## Literature boundary

Critical-angle measurements as independent bulk-refractive-index measurements and angular optimization of SPR detection predate CalSPR. V2.4 therefore does not claim invention of critical-angle refractometry, angular interrogation, or angle selection. Its bounded contribution is the equal-budget, equal-nuisance comparison and the link from the resulting observation design to a calibration-transfer risk threshold.

