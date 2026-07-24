# C-axis orientation distribution functions

## 1. Objective

Add manuscript-ready c-axis orientation distribution functions to the existing
comprehensive EBSD texture outputs. The new results must distinguish:

1. the acute polar angle between the alpha-Ti c-axis and the bar axial direction
   (AD); and
2. the complete c-axis distribution on the AD-positive projective hemisphere,
   including azimuth about AD.

The figures complement, but do not replace, the existing `{0001}` pole figures,
mean c-axis angle trend, and per-orientation
`c_axis_orientation_distribution.csv`.

## 2. Inputs and evidentiary hierarchy

Use the existing comprehensive EBSD texture data under
`results/mtex_ebsd_comprehensive/05_texture/`.

- Horizontal map direction and specimen x are AD.
- Specimen y is the in-plane transverse/radial direction (TD/RD).
- Specimen z is the scan normal (ND).
- The alpha-Ti c-axis is antipodal, so `c` and `-c` are equivalent.
- `raw_full` is the primary condition-level result.
- `denoised_full` is shown in separate figures as requested.
- `denoised_raw_common` is retained as a paired processing-sensitivity result
  and is not treated as an independent sample.
- Pixel/area weighting is the primary weighting.
- Area-weighted grain-mean orientations are a sensitivity result and do not
  appear as additional curves in the main figures.

No raw CTF or existing output is overwritten.

## 3. One-dimensional c-axis-to-AD distribution

### 3.1 Angle definition

For each c-axis unit vector `c` and AD unit vector `x`, calculate

```text
theta = arccos(abs(c dot x))
```

where `theta` is the acute c-axis-to-AD angle in the closed interval from
0 to 90 degrees.

### 3.2 Empirical probability density

The primary bin width is 2 degrees, with bin edges fixed from 0 to 90 degrees.
For each bin:

```text
probability = weighted count in bin / total valid weight
pdf_per_degree = probability / bin_width_deg
```

The probability over all bins must sum to one within numerical tolerance.
Sensitivity tables use 1, 2, and 5 degree bins.

### 3.3 Random-orientation normalization

A uniform distribution of antipodal axes on a hemisphere is not uniform in
`theta`; its probability is proportional to `sin(theta)`. The exact random-axis
probability of a bin with lower and upper bounds `theta_1` and `theta_2` is:

```text
random_probability = cos(theta_1) - cos(theta_2)
```

with angles evaluated in radians or with the corresponding degree-aware
functions.

The multiple-of-random-distribution value is:

```text
mrd = observed_probability / random_probability
```

A random axial distribution therefore has an expected MRD of 1 in every bin.
This bin-integrated definition avoids dividing a noisy point estimate by
`sin(theta)` near zero degrees.

### 3.4 One-dimensional figures

Create two separate figures:

```text
c_axis_ad_distribution_raw.png
c_axis_ad_distribution_denoised.png
```

Each figure contains:

- an upper panel showing the six probability-density curves in
  probability per degree;
- a lower panel showing the corresponding MRD curves;
- a horizontal MRD = 1 random baseline;
- the registered six-condition colour mapping and sample order;
- x-axis limits of 0–90 degrees;
- directly stated weighting, support, and 2 degree bin width.

Raw and denoised curves are not placed in the same main figure. Axis ranges are
kept identical between the two figures to permit visual comparison.

## 4. Complete spherical c-axis distribution

### 4.1 Hemisphere convention

Map each antipodal c-axis to the existing deterministic hemisphere convention:

1. AD-positive;
2. if the AD component is numerically zero, TD/RD-positive;
3. if both are numerically zero, ND-positive.

For the selected representative, calculate:

```text
theta_ad = arccos(c_x)
phi_about_ad = atan2(c_z, c_y)
```

where `theta_ad` spans 0–90 degrees and `phi_about_ad` spans -180–180 degrees.
The azimuth origin and positive direction are stated on every figure.

### 4.2 Spherical density

Estimate an antipodal c-axis density on a registered equal-area hemisphere grid.
The primary kernel halfwidth is 5 degrees, matching the current texture
workflow. The density is normalized so its hemisphere mean equals one and is
reported in MRD.

Sensitivity outputs compare kernel halfwidths of 5, 7.5, and 10 degrees. Grid
resolution is fine enough that changing the evaluation grid does not move a
stable density maximum by more than the registered tolerance.

The spherical density table records:

- sample and cold reduction;
- variant and support;
- weighting;
- kernel and grid parameters;
- `theta_ad` grid coordinate;
- `phi_about_ad` grid coordinate;
- equal-area cell weight;
- normalized MRD.

### 4.3 Spherical-distribution figures

Create two separate six-panel figures:

```text
c_axis_spherical_distribution_raw.png
c_axis_spherical_distribution_denoised.png
```

Each panel represents one deformation state in the registered order. Both
figures use:

- identical projection and orientation;
- AD, TD/RD, and ND annotations;
- one common MRD colour scale derived by a documented robust rule and registered
  for both variants;
- the same contour levels;
- the same sample labels and cold-reduction values;
- marked density maxima where numerically stable.

This output must not be presented as a full crystal ODF. It is a distribution of
one crystallographic axis and does not contain rotation about the c-axis.

## 5. Output tables and registration

Add the following artifacts to `05_texture`:

```text
c_axis_ad_distribution_function.csv
c_axis_spherical_distribution_function.csv
c_axis_distribution_parameters.csv
c_axis_ad_distribution_raw.png
c_axis_ad_distribution_denoised.png
c_axis_spherical_distribution_raw.png
c_axis_spherical_distribution_denoised.png
```

The one-dimensional CSV contains:

```text
sample
diameter_mm
cold_reduction_percent
variant
support
weighting
bin_width_deg
bin_lower_deg
bin_upper_deg
bin_center_deg
valid_source_count
valid_source_weight
observed_probability
pdf_per_degree
random_probability
mrd
```

The parameter CSV records coordinate definitions, hemisphere convention,
weighting, bin widths, kernel family, kernel halfwidths, spherical grid
resolution, colour-limit rule, primary support, and sensitivity supports.

The comprehensive output contract, artifact inventory, README, manifest
parameters, and corresponding tests are extended without changing existing
artifact definitions.

## 6. Implementation boundaries

The generator reads the existing per-orientation c-axis CSV and does not
re-import all CTF files when the required registered orientations and weights
are present.

For `denoised_raw_common`, pixel-weighted rows are paired by persistent source
ID and sample. Phase-changed or missing pairs are excluded and counted.
Grain-mean rows are not paired by grain ID across processing variants because
grain reconstruction may change their identity.

The one-dimensional distribution calculation and spherical density evaluation
are separate pure functions so they can be tested without plotting or file I/O.

## 7. Validation

### 7.1 Synthetic numerical tests

Test at least:

- all c-axes parallel to AD: all one-dimensional probability falls in the first
  bin and the spherical maximum is at AD;
- all c-axes normal to AD: all probability falls in the last included bin;
- antipodal inputs `c` and `-c`: identical one- and two-dimensional results;
- a deterministic random-axis quadrature: MRD approaches one across bins and on
  the spherical grid;
- a two-component mixture: the expected probability fractions and two
  spherical maxima are recovered;
- probability conservation;
- positive random-bin probabilities;
- hemisphere-mean spherical MRD equal to one;
- raw/common-mask source-ID pairing and exclusion counts.

### 7.2 Integration and visual checks

Verify:

- exactly six raw and six denoised states occur in each relevant output;
- all angles, probabilities, and MRD values are finite and within their defined
  ranges;
- each one-dimensional distribution sums to one;
- raw and denoised figure axes are identical;
- both spherical figures use the same registered colour scale;
- AD is oriented and labelled consistently;
- no labels, legends, or colour bars overlap;
- input and existing artifact hashes remain unchanged where registered;
- MATLAB `checkcode`, focused tests, output-contract tests, image inspection,
  and `git diff --check` pass.

## 8. Interpretation boundaries

The one-dimensional distribution can support statements about c-axis alignment
or redistribution in polar angle relative to AD. The spherical distribution can
separate polar-angle changes from azimuthal redistribution about AD.

These outputs do not directly demonstrate:

- continuous rotation of the same grains across deformation conditions;
- a specific active slip or twinning system;
- a temporal mechanism sequence;
- statistically significant process-wide differences from one field of view per
  condition.

A mean-angle change is interpreted together with peak position, peak width,
multimodality, and integrated component fractions. If the mean is nearly
constant while the distribution splits or component fractions exchange, the
manuscript describes texture redistribution rather than net c-axis rotation.

## 9. Completion criterion

The feature is complete when the two one-dimensional and two spherical figures,
their auditable CSV tables, parameters, contract registration, automated tests,
and visual checks are complete for all six raw and six denoised scans without
overwriting prior results.
