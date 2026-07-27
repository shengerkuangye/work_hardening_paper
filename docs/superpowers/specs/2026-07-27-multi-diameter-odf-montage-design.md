# Multi-diameter ODF montage design

## Purpose

Generate one manuscript-ready orientation distribution function (ODF) figure
that compares all six Gr4B23271 bar diameters. The figure must show complete
ODF sections rather than pole figures or the one-dimensional c-axis-AD angle
distribution.

## Selected design

The main figure uses a 3-row by 2-column montage. Each diameter occupies one
panel, and each panel contains the three Bunge Euler-space sections
`phi2 = 0 deg`, `30 deg`, and `60 deg`. The diameter sequence is 7, 6.48,
6.02, 5.6, 5.25, and 5 mm, following increasing cold reduction.

This layout is preferred over a 6-row by 1-column figure, which is too tall
for a manuscript page, and a 2-row by 3-column figure, which makes each
three-section ODF panel too narrow.

## Data scope and coordinate convention

- Input: the six registered raw CTF files in
  `data/ebsd_kpl_250221_7_df/scans`.
- Phase: indexed `Ti-Hex` pixels only.
- Primary statistical representation: pixel-weighted orientations.
- The denoised CTF files and area-weighted grain-mean ODF are excluded from
  this main figure; they remain sensitivity analyses in the comprehensive
  EBSD workflow.
- Specimen coordinates follow the registered project convention:
  X = AD, Y = TD/RD, and Z = ND.
- Sample labels include diameter and cold reduction:
  7 mm/0%, 6.48 mm/14.31%, 6.02 mm/26.04%, 5.6 mm/36.00%,
  5.25 mm/43.75%, and 5 mm/48.98%.

## ODF calculation

For each raw scan:

1. Load the CTF with the existing MTEX coordinate conversion.
2. Select valid `Ti-Hex` orientations.
3. Estimate the ODF with the De la Vallee Poussin kernel using a 5 deg
   halfwidth.
4. Normalize the ODF to a unit mean so intensity is reported in multiples of
   a random distribution (MRD).
5. Evaluate and render the three `phi2` sections with 5 deg angular
   resolution.

All six panels use one global color range from 0 to the maximum ODF intensity
evaluated across the six samples. A single shared colorbar is labelled
`ODF intensity (MRD)`. Sharing the scale is required for quantitative visual
comparison between diameters.

## Figure styling and outputs

- Layout: three rows by two columns, ordered left-to-right and top-to-bottom
  by decreasing diameter.
- Each panel title reports diameter and cold reduction.
- Section headings use `phi2 = 0 deg`, `30 deg`, and `60 deg`.
- The continuous color map must be perceptually ordered and must not use
  `jet` or a rainbow scale.
- Background is white and typography is consistent with the existing EBSD
  figures.
- The main raster output is exported at publication resolution, and a PDF
  companion is exported for manuscript assembly.

Outputs:

- `results/mtex_odf_diameter_montage/odf_diameter_montage.png`
- `results/mtex_odf_diameter_montage/odf_diameter_montage.pdf`
- `results/mtex_odf_diameter_montage/odf_diameter_summary.csv`

The summary table records sample name, diameter, cold reduction, input path,
valid Ti-Hex orientation count, kernel halfwidth, grid resolution, and the
maximum evaluated ODF intensity.

## Reproducibility and failure handling

The generator is a new focused MATLAB entry point under `tools/mtex`. It uses
`comprehensive_ebsd_catalog` and the registered loading helpers rather than
duplicating sample paths or coordinate transformations. It creates only
derived results and never modifies CTF inputs.

Generation stops with an explicit assertion if MTEX is unavailable, an input
file is missing, a scan contains no indexed `Ti-Hex` orientations, an ODF
normalization is invalid, or an expected output is empty.

## Verification

Automated checks must verify:

- exactly six raw catalog rows are selected in the registered diameter order;
- each scan contributes a finite, positive orientation count and ODF maximum;
- all six rendered panels use the same global MRD upper limit;
- the CSV contains six rows with the specified columns;
- PNG and PDF outputs exist and are non-empty.

The final PNG is visually inspected for readable labels, complete
`phi2` sections, correct diameter ordering, one shared colorbar, consistent
color limits, and absence of clipping or blank panels.
