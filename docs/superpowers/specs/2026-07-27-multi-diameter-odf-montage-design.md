# Multi-diameter ODF diagnostic design

## Purpose

Generate a seven-section orientation distribution function (ODF) diagnostic for each of the six Gr4B23271 bar diameters. This diagnostic supports subsequent evidence-based selection of manuscript ODF sections; it is not a pole-figure comparison or a one-dimensional c-axis--AD angle distribution.

## Selected design

The diagnostic is one 6-row by 7-column matrix. Rows are the registered bar diameters in this order: 7, 6.48, 6.02, 5.6, 5.25, and 5 mm. Columns are the Bunge Euler-space sections `phi2 = 0, 10, 20, 30, 40, 50, 60 deg`.

Each matrix cell is one `phi2` section; row labels report diameter and cold reduction, and column headings report `phi2`. The seven sections are the required diagnostic coverage of the alpha-Ti rotational fundamental zone. They must not be reduced to a sparse subset presented as the entire distribution.

No manuscript representative sections are fixed in this design. They shall be chosen only after the seven-section diagnostic has been generated and the per-section maxima have been inspected. The eventual manuscript caption shall describe them as **selected `phi2` sections** and state their actual angles.

## Data scope, symmetry, and coordinate convention

- Input: the six registered raw CTF files in `data/ebsd_kpl_250221_7_df/scans`.
- Phase: indexed `Ti-Hex` pixels only.
- Primary statistical representation: pixel-weighted orientations.
- The denoised CTF files and area-weighted grain-mean ODF are excluded from this diagnostic; they remain sensitivity analyses in the comprehensive EBSD workflow.
- Crystal symmetry: alpha-Ti is configured as `6/mmm`. ODF rotations are evaluated in the corresponding `622` rotational fundamental zone.
- Specimen symmetry: `SS = 1` for the bar specimens. Do not apply the plate-symmetry `222` assumption at any stage of ODF estimation, evaluation, or plotting.
- Specimen coordinates follow the registered project convention: X = AD, Y = TD/RD, and Z = ND.
- Sample labels include diameter and cold reduction: 7 mm/0%, 6.48 mm/14.31%, 6.02 mm/26.04%, 5.6 mm/36.00%, 5.25 mm/43.75%, and 5 mm/48.98%.

## ODF calculation and common comparison scale

For each raw scan:

1. Load the CTF with the existing MTEX coordinate conversion.
2. Select valid `Ti-Hex` orientations and retain pixel weighting.
3. Estimate the ODF with the same De la Vallee Poussin kernel for every sample, using a 5 deg halfwidth, `6/mmm` crystal symmetry, and `SS = 1`.
4. Normalize the ODF to a unit mean so intensity is reported in multiples of a random distribution (MRD).
5. Evaluate every `phi2 = 0:10:60 deg` section at 5 deg angular resolution.
6. Record the maximum MRD for every sample--section pair.

All 42 cells use one global color range from 0 to the largest evaluated MRD among all six samples and all seven sections. The same 5 deg kernel and 5 deg evaluation resolution apply to every cell. A single shared colorbar is labelled `ODF intensity (MRD)`.

## Figure styling and outputs

- Layout: six rows by seven columns, with rows ordered by decreasing diameter and columns ordered by increasing `phi2`.
- The continuous color map must be perceptually ordered and must not use `jet` or a rainbow scale.
- Background is white and typography is consistent with the existing EBSD figures.
- Export the full diagnostic at publication resolution as PNG and PDF.
- Do not export a manuscript representative-section figure until the full diagnostic and section-maxima CSV have been inspected and the selected `phi2` sections are documented.

Outputs:

- `results/mtex_odf_diameter_montage/odf_diameter_full_sections.png`
- `results/mtex_odf_diameter_montage/odf_diameter_full_sections.pdf`
- `results/mtex_odf_diameter_montage/odf_diameter_summary.csv`
- `results/mtex_odf_diameter_montage/odf_diameter_section_summary.csv`

`odf_diameter_summary.csv` has exactly six rows, one per diameter. It records sample name, diameter, cold reduction, input path, valid Ti-Hex orientation count, kernel halfwidth, grid resolution, the largest of that sample's seven section maxima, and the shared global color limit.

`odf_diameter_section_summary.csv` has exactly 42 rows (six samples times seven `phi2` sections). It records the preceding sample metadata together with `phi2_deg`, the section maximum MRD, and the shared global color limit. This table is the explicit selection record used to identify candidate manuscript sections after visual inspection.

## Reproducibility and failure handling

The generator is a new focused MATLAB entry point under `tools/mtex`. It uses `comprehensive_ebsd_catalog` and the registered loading helpers rather than duplicating sample paths or coordinate transformations. It creates only derived results and never modifies CTF inputs.

Generation stops with an explicit assertion if MTEX is unavailable, an input file is missing, a scan contains no indexed `Ti-Hex` orientations, an ODF normalization is invalid, the configured crystal symmetry is not `6/mmm`, the specimen symmetry is not `SS = 1`, an expected section is missing, or an expected output is empty.

## Verification

Automated checks must verify:

- exactly six raw catalog rows are selected in the registered diameter order;
- `6/mmm` crystal symmetry, its `622` rotational fundamental zone, and `SS = 1` are used; no `222` specimen symmetry is accepted;
- each scan contributes a finite, positive orientation count;
- all 42 sample--section pairs have a finite, positive section maximum;
- every cell uses the same global MRD upper limit, equal to the maximum across all 42 recorded section maxima;
- the six-row and 42-row CSV files have the specified columns;
- PNG and PDF outputs exist and are non-empty.

The final PNG is visually inspected for the six ordered rows, seven `phi2` columns, readable labels, a shared MRD colorbar, identical color limits, and absence of clipping or blank cells. Only after that inspection may a manuscript figure be designed from selected `phi2` sections.
