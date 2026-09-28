# CalSPR V2.4 reproducibility archive

Companion computational materials for "Information-source and risk boundaries for
SPR calibration transfer". This archive contains planar optical and count models,
inversion, Fisher-information design, predeclared configurations, deterministic tests,
material-source data, frozen numerical outputs and figure source tables. These are
computational results, not measurements on a manufactured instrument.

## Contents and preservation

scripts/calspr contains the complete numerical implementation, including inversion.py
and v24_design.py. data/raw/materials contains all three original CC0 optical-constant
tables. configs retains frozen settings; outputs retains final results and required
source stages. Figures are regenerable from the supplied CSV files. No manuscript,
author email list, credentials, original execution logs or binary dependency caches
are included. Source numeric files and material tables are byte-preserved; only the
documented root/config path expressions in release helper scripts were adapted.

SHA256SUMS.txt records archive content. verification/SOURCE_FILE_HASHES.json records
the original/copied file hashes, and PORTABILITY_CHANGES.json records path changes.
The three required core files and material files can be used after moving this root.

## Environment

The validated research host uses Linux x86_64, Python 3.12.3 and CPU computation.
Install dependencies in an isolated environment using requirements.txt; see
verification/DEPENDENCY_VERSIONS.json for actual versions verified during packaging.
No CUDA is used. Set OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1 and
NUMEXPR_NUM_THREADS=1. Predeclared parallel confirmations use at most four workers.

The plotting scripts use DejaVu Sans. Pillow helpers currently refer to the common
Linux DejaVu system-font paths; font portability is a documented plotting dependency.

## Quick deterministic verification

From this directory, run:

```sh
python3 tests/run_tests.py
```

This overwrites only outputs/test_summary.json in your working copy. For complete
run logging use `bash scripts/run_logged.sh tests configs/p0_p2.json 0`. The test
algorithm is deterministic (seed not applicable; 0 is only PYTHONHASHSEED).
Tests cover interfaces, energy balance, independent optical solvers, limiting cases
and the SPR feature. The historical test script records a tests=11 field but emits
12 PASS checks; validation logs record the actual count.

## Existing result map

- Figure 1-5 source tables: figures/source_data and outputs/boundary_scan,
  outputs/boundary_confirmation, outputs/target_recalibration_* and integrity_audit.
- Figure 6 information design: outputs/v24/v24_scan_20260914_e4_final and
  outputs/v24/v24_design_confirmation_20260914_c1 (100 replicates per design/state).
- Figure 7 transfer-scale confirmation: outputs/v24/v24_confirmation_20260914_c3_combined_n200
  (200 paired replicates per condition; 25,200 strategy records). The c1 base and
  c2 extension are retained to trace the predeclared extension.

Exact seeds are in configs and each output summary: the V2.4 scan uses 2026091401,
paired scale confirmation 2026091402, and design confirmation 2026091403. Earlier
boundary and target-recalibration runs have distinct seeds; do not merge their designs.

## Rerunning analyses

Frozen run directories are present, and analysis scripts refuse to overwrite them.
For independent reruns, create a separate working copy and move that copy's outputs
to an archival sibling directory first; never delete this archive's evidence.
Keep required upstream input results at their recorded paths when rerunning downstream
stages. Use unique run IDs when exploring changed designs, and label them as new runs.

```sh
bash scripts/run_logged.sh v24_scan configs/v24_scan_split_final.json 2026091401
bash scripts/run_logged.sh v24_confirmation configs/v24_confirmation_frozen_fast_v1.json 2026091402
bash scripts/run_logged.sh v24_design_confirmation configs/v24_design_confirmation_frozen_v1.json 2026091403
python3 scripts/release/run_v24_confirmation_extension.py configs/v24_confirmation_extension_v1.json
python3 scripts/release/merge_v24_confirmation_n200.py
```

The last two commands require the base run and frozen trigger/source hashes. They
must be executed with equivalent command/time/version/seed/input/output/exit logging.
Full Monte Carlo was not rerun during archival preparation. Existing optimizer
checks/results are retained; packaging verification is not a replacement Monte Carlo.

For plotting, `python3 scripts/release/make_v22_main_figures_final_author_cleanup.py
--help` shows project/output-root arguments. The V2.4 Pillow helpers accept the
corresponding configs/v24_figures_pillow_v3_final.json or
configs/v24_figure7_final_author_cleanup_20260917.json. Change output_dir and
source_data_dir in a copied figure config to an unused location before plotting.
Final publication layout may differ from the archived research rendering.

## Release status

This archive is prepared for public deposit under the author-selected MIT
license for original CalSPR code. Third-party plotting helpers remain Apache-2.0
and optical-constant tables retain CC0. See LICENSE_STATUS.md and THIRD_PARTY_NOTICES.md.
Source repository: https://github.com/mzw423502/calspr-calibration-transfer

Software release version: 2.4.1 (scientific study version V2.4). The GitHub
repository has been created; this package supplies its release content. An
archival DOI has not yet been assigned. CITATION.cff and .zenodo.json provide
matching citation and deposit metadata without a fabricated identifier.

## Measurement revision figure numbering (2026-09-28)

Historical to submission figure mapping is 6->1, 2->2, 7->3, 1->4, 3->5, 4->6,
5->7 (see FIGURE_NUMBER_MAPPING.json). Historic script function names are preserved.
scripts/release/make_measurement_figures_1_to_5.py contains the supplied 2026-09-21
layout fix for historical Figure 5 (now Figure 7). From the archive root, use
the actual positional arguments and a fresh output directory:

```sh
python3 scripts/release/make_measurement_figures_1_to_5.py . figures/measurement_rerender
```

Then rename/copy exports using FIGURE_NUMBER_MAPPING.json; never rewrite source values. Final 2026-09-28 Figure 1/3 display repairs, original raster inputs,
frozen CSVs and verified rendering/audit artifacts are included under
figures/revision_20260928; use scripts/figures/fix_frozen_figures.py. This Pillow
rendering step was verified separately on the research host.

The archive preparation host cannot currently import the optional vendored
Matplotlib binary (shared-object load restriction); numerical imports and all 12
deterministic optical checks passed. Figure code compiled, but complete Matplotlib
rendering was not claimed or performed during this package check. Existing final
tables and visual files in the manuscript package remain the scientific source.
