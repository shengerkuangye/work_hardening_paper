# OIM-style square ODF peak-matrix design

## Purpose

Create one publication-oriented ODF matrix that reproduces the square-section
presentation of the supplied OIM Analysis figure while preserving the
project's unsymmetrized quantitative ODF convention. The matrix compares all
six Gr4B23271 diameters at seven Bunge Euler-space sections and identifies the
strongest orientation visible in every panel.

The new result supplements, and does not overwrite or replace, the existing
full-`phi1` ODF diagnostic. Because the square display is a restricted window
of an `SS = 1` ODF, its annotated peak is explicitly reported as a
display-region peak rather than the global maximum of the complete section.

## Selected interpretation of the OIM workflow

The ODF calculation retains specimen symmetry `SS = 1`. No orthorhombic
specimen symmetry (`222`/`mmm`) is imposed. This follows the stated OIM
workflow, in which a square angular range is selected for the ODF output but
no specimen-symmetry operation is specified.

- The hexagonal crystal symmetry restricts the unique `phi2` interval to
  `0--60 deg`.
- With `SS = 1`, the complete `phi1` interval remains `0--360 deg`.
- The plotted square `phi1 = 0--90 deg`, `Phi = 0--90 deg` is therefore a
  display window, not a complete specimen-symmetry fundamental region.

This distinction must appear in the generated parameter record and in any
manuscript caption that uses the figure.

## Input data and ODF model

- Input: the six registered raw CTF files in
  `data/ebsd_kpl_250221_7_df/scans`.
- Sample order: 7, 6.48, 6.02, 5.60, 5.25, and 5.00 mm.
- Phase: indexed `Ti-Hex` pixels only.
- Statistical representation: pixel-weighted orientations.
- Crystal symmetry: imported `6/mmm`, with proper rotational group `622`.
- Specimen symmetry: imported `SS = 1`.
- Specimen coordinates: X = AD, Y = TD/RD, and Z = ND.
- ODF kernel: De la Vallee Poussin, 5 deg halfwidth.
- Density normalization: unit mean, reported as multiples of a random
  distribution (MRD).

The generator shall reuse the registered catalog and loading helpers. It shall
not modify the raw CTF data or alter the existing ODF generators and results.

## Section evaluation and peak definitions

The seven plotted columns are fixed at
`phi2 = 0, 10, 20, 30, 40, 50, and 60 deg`.

For each sample and `phi2` section, evaluate the harmonic ODF on a 1 deg Bunge
Euler grid using two domains:

1. **Square display domain:** `phi1 = 0--90 deg` and `Phi = 0--90 deg`,
   including both endpoints. Its maximum is the value marked in the panel.
2. **Complete `SS = 1` section:** `phi1 = 0--359 deg` and
   `Phi = 0--90 deg`. Its maximum is recorded as a quantitative audit value
   but is not marked if it lies outside the displayed square.

Record the MRD value and Bunge coordinates `(phi1, Phi, phi2)` for both
maxima. Also record whether the complete-section maximum lies in the square
display domain and the MRD difference between the two maxima.

Grid ties are resolved deterministically by the evaluation order: lowest
`Phi`, then lowest `phi1`. The recorded coordinates are grid estimates with
1 deg angular precision; they are not presented as exact crystallographic
components or evidence of a particular slip system.

## Orientation description

Each display-region peak receives a stable identifier `P01` through `P42`,
ordered first by decreasing diameter and then by increasing `phi2`.

The panel contains:

- a high-contrast peak marker at the evaluated `(phi1, Phi)` location;
- the peak identifier;
- the display-region peak intensity in MRD.

The companion peak table maps every identifier to diameter, cold reduction,
`phi2`, display-region peak MRD, and its Bunge Euler angles. It also reports
the complete-section peak and audit fields described above.

For an interpretable but non-speculative texture note, calculate the acute
angles between the crystal c-axis `[0001]` at the display-region peak and the
three registered specimen axes AD, TD/RD, and ND. A generated text field names
only the nearest specimen axis, for example
`c-axis closest to TD/RD (18.4 deg)`. No basal, prismatic, pyramidal, fibre, or
slip-system label is assigned automatically.

## Matrix layout and styling

- Layout: six rows by seven columns.
- Rows: diameters in decreasing order, labelled with diameter and registered
  cold-reduction percentage.
- Columns: increasing `phi2`, labelled at the top.
- Every panel is geometrically square and spans
  `phi1 = 0--90 deg`, `Phi = 0--90 deg`.
- Outer-axis labels identify horizontal `phi1` and vertical `Phi`; sparse
  `0`, `45`, and `90 deg` ticks keep the 42-panel matrix legible.
- All panels use identical filled-contour levels and one shared colorbar
  labelled `ODF intensity (MRD)`.
- The shared upper color limit is the largest evaluated square-domain value
  across the 42 plotted panels. Thus no plotted peak is clipped and visual
  comparison between diameters and sections is quantitative.
- Plotting-only negative interpolation artefacts are clipped to `0 MRD`;
  ODF models and tabulated values remain unchanged.
- Use a perceptually ordered colormap rather than `jet` or an OIM-style
  rainbow scale. Contour boundaries and peak markers remain distinguishable
  in print.
- Background, typography, margins, and panel geometry are uniform.

The main title and parameter record state `SS = 1` and
`display window: phi1 = 0--90 deg`. A concise figure note states that marked
peaks are maxima within the displayed square region.

## Outputs

Create a separate derived-results directory:

`results/mtex_odf_square_peak_matrix/`

The generator produces:

- `odf_square_peak_matrix.png`
- `odf_square_peak_matrix.pdf`
- `odf_square_peak_summary.csv`
- `odf_square_peak_parameters.csv`

The PNG and image-content PDF are exported at approximately 600 dpi.

`odf_square_peak_summary.csv` contains exactly 42 rows and includes:

- peak identifier;
- sample, diameter, cold reduction, and input path;
- `phi2`;
- display-domain peak MRD and Bunge `phi1`, `Phi`, `phi2`;
- c-axis angles to AD, TD/RD, and ND;
- nearest-axis texture note;
- complete-section peak MRD and Bunge `phi1`, `Phi`, `phi2`;
- whether the complete-section maximum is inside the display domain;
- complete-minus-display peak MRD;
- shared color limit.

`odf_square_peak_parameters.csv` is a one-row reproducibility record containing
the data variant, phase, weighting, crystal and specimen symmetries, coordinate
convention, calculation method, kernel type and halfwidth, grid resolution,
section angles, square and complete evaluation domains, normalization,
tie-breaking rule, and color-limit definition.

## Failure handling

Generation stops with an explicit error if:

- MTEX or a required project helper is unavailable;
- any registered input file is absent;
- the six samples are missing, duplicated, or out of order;
- a scan has no indexed `Ti-Hex` orientations;
- imported or harmonic symmetries differ from `6/mmm`, `622`, and `SS = 1`;
- an ODF cannot be normalized to a finite positive unit mean;
- a section grid contains nonfinite values;
- any peak coordinate lies outside its declared domain;
- the expected 42 peak records are not produced;
- a required output is empty.

## Verification

Automated tests shall verify:

- six raw scans in the registered order and seven sections exactly
  `0:10:60 deg`;
- no `222`/`mmm` specimen symmetry is introduced;
- the square and complete domains contain the expected grid coordinates;
- all 42 display and complete-section peaks are finite and positive;
- display peaks lie within `0--90 deg` for both plotted Euler axes;
- complete-section peaks lie within `0--359 deg` in `phi1` and
  `0--90 deg` in `Phi`;
- every complete-section peak is at least as large as its corresponding
  display-region peak, within numerical tolerance;
- peak identifiers are unique and ordered `P01` through `P42`;
- c-axis angles are finite and the nearest-axis note agrees with their
  minimum;
- the shared color limit equals the largest display-domain peak;
- the peak and parameter CSV schemas and row counts are exact;
- PNG and PDF outputs exist, are non-empty, and retain approximately 600 dpi
  image content.

Final visual inspection shall confirm:

- a complete 6-by-7 matrix with square panels;
- row and column order;
- readable peak markers, identifiers, MRD labels, and outer Euler-axis labels;
- one common color scale with no clipped plotted peaks;
- no blank panels, contour gaps, label collisions, or cropped margins.

The manuscript caption must call the sections
`selected phi2 sections of the SS = 1 ODF` and the marked values
`maxima within the displayed phi1 = 0--90 deg and Phi = 0--90 deg region`.
