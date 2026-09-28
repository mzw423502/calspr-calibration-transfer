# CalSPR V2.4 numeric and source index

| Claim or artifact | Value/status | Authoritative source |
|---|---|---|
| Critical angle | 61.6206773 degrees | `outputs/v24/v24_scan_20260914_e4_final` |
| Frozen windows | C 60.7-62.6 (20); S 71-75 (41); O 140; F 201 | scan output and `docs/v24/V24_MEASUREMENT_AND_WINDOWS.md` |
| Projected shares | C 0.252181%; S 77.094623%; O 22.653195% | `outputs/v24/v24_final_audit_20260914_v2/V24_CLAIM_VALUES.json` |
| Equal-budget local SD | F 1.6712e-5; C 1.4917e-4; S/CS 8.7493e-6; P 2.7328e-5 RIU | same claim-values file |
| A ranking | S first by RMSE in 4/4 states | `V24_A_DESIGN_RANKING.csv` in final audit |
| Threshold range | 0.026721-0.429103 nm | final audit and `figures/v24_final/source_data/V24_THRESHOLD_PREDICTIONS.csv` |
| 50 nm, 1e7 finite threshold | 0.031024585 nm | scan output |
| n=200 fit count | 25,200 | combined B replicate file |
| Risk direction | zero 9/9 T0 mean; half 16/18 T0 mean; twice 18/18 T1 mean; one 17/18 unresolved | `V24_RISK_DIRECTION_AUDIT.csv` |
| Residual precision triggers | 14 | `docs/v24/V24_RESIDUAL_PRECISION_STATES.csv` |
| Technical final audit | PASS 23/23 | `outputs/v24/v24_final_audit_20260914_v2/summary.json` |
| Literature boundary | no priority claim | `docs/v24/V24_LITERATURE_POSITIONING.md` |

Complete file hashes are supplied by the release manifest; this table is a claim-to-source map, not a substitute for it.
