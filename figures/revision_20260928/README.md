# Verified 2026-09-28 display repair

The inputs, exact Pillow script and already verified PNG/PDF/audit output are
included. Historical Figure 6 becomes submission Figure 1; historical Figure 7
becomes submission Figure 3. Other historical panels remain pixel-preserved except
the repaired 6c/6d and 7c and the documented title text. Frozen CSVs supply the
newly drawn axes and marks; no model fitting or Monte Carlo is rerun.

From the archive root (Pillow and Linux DejaVu Sans required):

    python3 scripts/figures/fix_frozen_figures.py --input-dir figures/revision_20260928/inputs --output-dir figures/revision_20260928/new_render

Use a fresh output directory. verified_render/plot_data_audit.json records exact
input rows; text_bbox_audit.json records layout boxes; manifest.json records hashes.
The supplied research-host log stays in the private project; public inputs and
outputs are sufficient to rerun this deterministic rendering step.
