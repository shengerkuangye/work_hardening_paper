# Selected-phi2 ODF figure design

## Purpose

Create a presentation- and manuscript-oriented ODF figure for the six
Gr4B23271 bar diameters following the supervisor's requested interpretation:
quantify texture strength with the ODF maximum, identify texture components
from selected `phi2` sections, and retain peak positions for subsequent
crystal-rotation and slip-system analysis.

The existing seven-section diagnostic remains unchanged. The new figure is a
selected-section view and must not be described as the complete ODF.

## Data and ODF model

- Use the six registered raw CTF scans in the order 7, 6.48, 6.02, 5.60,
  5.25, and 5.00 mm.
- Select indexed `Ti-Hex` pixels and retain pixel weighting.
- Reuse the existing ODF settings: imported `6/mmm` crystal symmetry,
  rotational group `622`, specimen symmetry `SS = 1`, and a 5 deg
  De la Vallee Poussin kernel.
- Normalize the radial-basis and harmonic ODFs to unit mean with the guarded
  normalization helper.
- Do not apply specimen symmetry `222`.

## Quantitative texture intensity

For each sample, calculate the ODF maximum with MTEX
`max(odf,"resolution",1*degree)`. This value is the formal ODF maximum used
in the row annotation and summary table; it replaces the seven-section
5 deg grid maximum as the quantitative texture-intensity measure.

Record the orientation at the maximum in Bunge Euler angles
`phi1_max_deg`, `Phi_max_deg`, and `phi2_max_deg`. These coordinates support
later evaluation of crystal rotation but are not interpreted as proof of a
specific slip system.

All 18 selected panels use one common color range from zero to the largest
formal ODF maximum among the six samples. The colorbar is labelled
`ODF intensity (MRD)`.

## Selected sections and layout

- Layout: six rows by three columns.
- Columns: `phi2 = 0, 30, and 60 deg`.
- Rows: decreasing diameter in the registered order.
- Every panel spans the full Bunge ranges `phi1 = 0--360 deg` and
  `Phi = 0--90 deg`.
- Each left-side row label contains:
  - diameter;
  - registered cold-reduction percentage;
  - `Max ODF = x.xx MRD`.
- Column headings state the selected `phi2` angle.
- A single shared parula colorbar is placed at the far right.
- Background is white, contour lines are black, and typography is consistent
  across all panels.

Nonphysical negative values on the MTEX contour plotting grid are clipped to
`0 MRD` for display only. The ODF model, formal maximum, section peak values,
and CSV records remain unchanged.

## Peak-position records

In addition to the six-row ODF-maximum summary, evaluate each selected section
on a 1 deg `phi1`--`Phi` grid and record:

- sample and diameter;
- cold-reduction percentage;
- selected `phi2`;
- section peak MRD;
- `phi1` and `Phi` coordinates of the section peak;
- the common global color limit.

The 18-row selected-section table provides a reproducible basis for comparing
peak movement. Texture-component assignment and slip-system interpretation
remain a later analysis step.

## Outputs

Create a new derived-results directory without overwriting the seven-section
diagnostic:

- `results/mtex_odf_selected_phi2/odf_selected_phi2_sections.png`
- `results/mtex_odf_selected_phi2/odf_selected_phi2_sections.pdf`
- `results/mtex_odf_selected_phi2/odf_selected_phi2_summary.csv`
- `results/mtex_odf_selected_phi2/odf_selected_phi2_peak_positions.csv`

The PNG and PDF are exported at approximately 600 dpi. The PDF uses explicit
600 dpi image content to avoid MATLAB's default downsampling of composed
raster rows.

## Verification

Automated checks verify:

- six samples in the registered order;
- selected sections exactly `[0,30,60] deg` for every sample;
- six formal ODF maxima are finite, positive, and computed by MTEX at
  1 deg resolution;
- 18 section peaks are finite, positive, and have valid Euler coordinates;
- actual imported and harmonic ODF symmetries are `6/mmm`, `622`, and
  `SS = 1`;
- one common color limit equals the largest formal ODF maximum;
- contour plotting data are finite and nonnegative after display-only
  clipping;
- the two CSV schemas and row counts are exact;
- PNG and PDF are non-empty, approximately 600 dpi, and horizontally
  oriented.

Visual inspection verifies:

- a complete 6-by-3 matrix;
- column order `0, 30, 60 deg`;
- readable left-side diameter, reduction, and maximum-ODF annotations;
- no internal white contour gaps;
- one shared colorbar;
- no clipping or blank panels.

The figure caption must call these **selected `phi2` sections** and list the
three angles.
