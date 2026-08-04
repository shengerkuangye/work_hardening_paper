# Six-state alpha-Ti c-axis pole-figure montage design

## Purpose

Create one publication-ready figure that compares the alpha-Ti `{0001}`
c-axis pole-figure distributions for all six cold-reduction states. The figure
is intended to evaluate population-level c-axis redistribution; it must not be
presented as direct tracking of individual-grain rotation.

## Scientific data contract

- Input: the six raw EBSD CTF scans registered by
  `comprehensive_ebsd_catalog`.
- Phase: valid `Ti-Hex` indexed points only.
- Weighting: one equal weight per EBSD sampling point, equivalent to the
  existing pixel-weighted ODF/PF convention.
- ODF kernel: de la Vallee Poussin, 5 degree halfwidth.
- Pole: `{0001}` with hexagonal crystal symmetry and antipodal treatment.
- Projection: equal-area upper-hemisphere pole figure.
- Specimen coordinates: AD = sample X and horizontal, TD/RD = sample Y, and
  ND = sample Z.
- All six panels must share one MRD range calculated from the maximum across
  the six pole densities. Individual-panel auto-scaling is prohibited.

## Figure design

- Layout: two rows by three columns, ordered by cold reduction:
  `0`, `14.31`, `26.04`, `36.00`, `43.75`, and `48.98%`.
- Each panel title contains diameter and cold reduction.
- Direction annotations identify AD, TD/RD, and ND.
- A single shared colorbar is labelled `MRD`.
- The title identifies the figure as pixel-weighted alpha-Ti `{0001}` c-axis
  pole figures.
- Use a perceptually ordered continuous colormap and contour boundaries that
  remain readable in print. Do not use a separate legend or color scale per
  panel.

## Outputs

Write new files without replacing the six existing single-panel figures:

- `results/mtex_c_axis_pf/c_axis_pf_six_state_montage.png`
- `results/mtex_c_axis_pf/c_axis_pf_six_state_montage.pdf`
- `results/mtex_c_axis_pf/c_axis_pf_six_state_montage_metadata.csv`

The metadata records the six source scans, pixel counts, kernel halfwidth,
projection, pole, specimen-axis convention, each panel maximum, and the shared
MRD maximum.

## Implementation boundary

Add a focused MTEX generator and its test. Reuse the existing catalog, EBSD
loader, and texture conventions. Do not alter raw CTF files or overwrite the
existing individual PF outputs.

## Verification

- Automated test confirms six ordered panels, six valid raw inputs, positive
  Ti-Hex pixel counts, one shared finite MRD maximum, and all three output
  files.
- PNG is at least 600 dpi and has a landscape 2x3 layout.
- PDF is non-empty and begins with a valid PDF header.
- MATLAB Code Analyzer reports no issues in the new generator and test.
- Visual inspection confirms no clipping, correct direction labels, readable
  panel titles, one colorbar, and identical color limits in all panels.

## Interpretation boundary

The montage may support statements about c-axis distribution, broadening,
peak relocation, or component redistribution. Without same-area before/after
grain tracking, it cannot prove the rotation trajectory of an individual
grain.
