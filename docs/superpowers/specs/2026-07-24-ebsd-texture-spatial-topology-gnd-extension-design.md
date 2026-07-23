# EBSD texture, spatial-localization, LAGB-topology, and minimum-energy GND extension

## 1. Purpose and scope

This extension adds four evidence-focused modules to the existing comprehensive
EBSD workflow:

1. crystallographic texture-component and peak evolution;
2. spatial localization of KAM and GROD relative to boundary networks;
3. topology of the 2–5 degree and 5–15 degree low-angle boundary networks;
4. a conventional two-dimensional EBSD minimum-energy GND estimate constrained
   by observable curvature.

The extension is intended to determine what microstructural changes are directly
supported by the existing EBSD scans, without requiring the data to conform to a
preassigned sequence of “c-axis rotation first, dislocation multiplication
later”. That sequence remains a scientifically useful hypothesis. It will be
accepted, modified, or rejected according to the combined texture, spatial
orientation-gradient, boundary-network, and mechanical evidence.

This work extends, but does not replace, modules `00_audit` through
`09_sensitivity` under `results/mtex_ebsd_comprehensive/`. Module
`09_sensitivity` currently exists on disk but is not registered by the output
contract or main runner. Registration of Module 09 is a prerequisite integration
fix, not a new scientific analysis.

## 2. Inputs, geometry, and evidentiary hierarchy

The inputs remain the six original CTF scans and their six paired denoised CTF
files registered by the current comprehensive workflow.

- The horizontal map direction is the bar axial direction (AD), equal to the
  MTEX specimen x direction.
- The vertical map direction is the in-plane transverse/radial direction
  (TD/RD), equal to specimen y.
- The scan normal is specimen z.
- The native grid is 600 by 600 pixels at 0.5 micrometre spacing.
- Original CTF data are the primary quantitative evidence.
- Denoised data are paired sensitivity data and are not independent replicates.
- Raw input files are read-only and their SHA-256 hashes must remain unchanged.
- Calculations use indexed Ti-Hex pixels unless an output explicitly states a
  different population.

The primary interpretation is based on trends that are:

1. present in the original data;
2. not reversed by the denoised-data sensitivity result;
3. stable under the declared parameter sensitivity tests; and
4. not attributable to indexed-fraction changes, map edges, or isolated
   low-quality pixels.

Because denoising fills pixels that are unindexed in the original scans, every
raw-versus-denoised analysis uses three explicit support labels:

- `raw_full`: all indexed Ti-Hex pixels in the original scan, used for the
  primary condition trend;
- `denoised_raw_common`: denoised orientations restricted to coordinates that
  are Ti-Hex in both the original and denoised scans, used for the primary
  paired processing-sensitivity comparison;
- `denoised_full`: all indexed Ti-Hex pixels in the denoised scan, used
  separately to quantify the effect of filling/recovery.

Full-map raw and full-map denoised differences are not described as processing
robustness unless the common-support comparison agrees. The support audit
reports the indexed, common, filled, removed, and phase-changed populations.

## 3. Integration architecture

Create four additive output modules:

```text
results/mtex_ebsd_comprehensive/
  10_texture_components/
  11_spatial_localization/
  12_lagb_topology/
  13_gnd_lower_bound/
```

The modules must use the current sample registry, import convention, grain
reconstruction, manifest logic, figure style, and output-contract mechanism.
Shared definitions are centralized; module generators must not silently redefine
AD, sample order, reduction, reconstruction thresholds, or raw/denoised roles.

The existing output contract and verifier depend on directory/artifact field
order. Directories 09–13 and their artifact fields must therefore be appended in
the same explicit order and covered by a contract test. Per-state map artifacts
use a declared 12-file pattern rather than a single placeholder filename.

Each module writes:

- analysis parameters;
- a condition-level summary;
- any distribution or component-level table required to reproduce a plotted
  trend;
- individual or montage figures;
- raw-versus-denoised comparison fields;
- a machine-checkable coverage or validity count.

The existing `analysis_manifest.csv` is extended with the four module runs,
parameter-set identifiers, artifact paths, input hashes, and completion status.
It records input hashes and registered parameters; it is not described as an
output-artifact hash manifest unless output hashes are deliberately added.

The fresh-run order is:

1. generate the existing 00–06 modules;
2. generate and register 09 sensitivity;
3. generate 10 texture components;
4. generate 11 spatial localization;
5. generate 12 LAGB topology;
6. generate 13 minimum-energy GND estimates;
7. generate 07 raw/denoised comparisons;
8. generate 08 tensile integration;
9. write the manifest and README, then verify the complete contract.

This order allows only a concise, preregistered subset of the new summary
metrics to enter generic raw/denoised and tensile integration. Component,
distance-bin, graph-component, histogram, and sensitivity long tables do not
enter the generic comparison logic because their numeric identifier columns can
be misclassified as measurements.

## 4. Module 10: texture-component and peak evolution

### 4.1 Scientific question

The module distinguishes among:

- continuous displacement of a texture maximum;
- redistribution between two or more texture components;
- peak broadening;
- peak splitting or merging;
- a change in texture intensity with little change in peak location.

These alternatives cannot be separated reliably by a single mean c-axis angle.

### 4.2 Orientation population and ODF

The primary analysis uses pixel/area weighting because all native pixels have the
same mapped area. The existing area-weighted grain-mean ODF is retained as a
sensitivity analysis of the influence of intragranular orientation spread and
must be labelled explicitly.

The primary ODF retains the current common settings:

- SO(3) de la Vallée Poussin kernel;
- 5 degree halfwidth;
- 2.5 degree common SO(3) peak-evaluation grid.

The same crystallographic symmetry, specimen reference frame, grid, and
normalization are used for all twelve scans.

The 5 degree evaluation grid already used for global indices and plotting is
retained for backward compatibility. Peak location is evaluated at 2.5 degrees
because a 5 degree grid is too coarse for quantitative peak tracking. Evaluation
grids of 1, 2.5, and 5 degrees and ODF halfwidths of 5, 7.5, and 10 degrees form
the registered peak-stability sensitivity set.

### 4.3 Peak detection

For each scan and weighting:

1. evaluate the normalized ODF on the registered common SO(3) grid;
2. identify symmetry-aware local maxima;
3. merge maxima that represent the same symmetry-equivalent orientation;
4. retain peaks above a documented intensity and prominence rule;
5. report non-retention rather than forcing every condition to have the same
   number of peaks.

For each retained peak, report:

- component identifier within the scan;
- representative Euler angles and orientation object;
- peak MRD;
- local prominence;
- peak spread defined separately from the ODF kernel halfwidth, using a
  documented local weighted RMS misorientation and/or half-maximum support
  volume;
- integrated orientation fraction within 10 and 15 degrees;
- c-axis angles to AD, TD/RD, and ND;
- azimuth about AD under the existing hemisphere convention.

Euler angles are presentation coordinates only. Peak matching and angular
distances use symmetry-reduced crystallographic misorientation. The output
separates `{0001}` pole-figure maxima, which describe c-axis directions, from
full ODF maxima, which also contain rotation about the c-axis.

If 10 or 15 degree peak neighbourhoods overlap, the output reports both the
non-exclusive neighbourhood integral and an exclusive nearest-peak assignment.
Non-exclusive integrals are never summed and called total component fraction.

### 4.4 Peak tracking between deformation conditions

Candidate matches are constructed between adjacent reduction states using
symmetry-reduced orientation distance and similarity of peak support. Global
assignment uses a declared maximum angular gate. A one-to-one assignment is used
only when unambiguous. The output explicitly
supports:

- continuation;
- birth;
- disappearance;
- split;
- merge;
- ambiguous match.

No component is forced to continue across a crystallographically implausible
angular jump. The tracking table records the angular distance and match status.
Because the six conditions are not in-situ observations of the same grains, a
tracked peak is described as a condition-to-condition texture-peak trajectory,
not direct observation of continuous rotation of the same grain population.

Literature component names may be attached only when:

1. a primary source gives an explicit ideal orientation or fibre definition for
   commercially pure alpha titanium or a directly comparable HCP system; and
2. the measured peak falls within a declared tolerance.

Otherwise, components retain neutral data-driven identifiers.

### 4.5 Required outputs

```text
10_texture_components/
  texture_component_parameters.csv
  texture_component_peaks.csv
  texture_component_tracks.csv
  texture_component_sensitivity.csv
  texture_component_trends.png
  texture_component_peak_maps.png
```

`texture_component_peaks.csv` includes sample metadata, variant, weighting,
peak identity, peak orientation, MRD, prominence, width, 10/15 degree fractions,
and c-axis descriptors.

`texture_component_tracks.csv` includes source and destination condition,
component identifiers, symmetry-reduced angular distance, match status, and
changes in MRD, width, and integrated fraction.

### 4.6 Interpretation boundary

A systematic change in the `{0001}` peak direction supports net c-axis
reorientation. A full ODF-peak displacement without a corresponding `{0001}`
change may instead reflect rotation about the c-axis. Exchange of integrated
fraction between stationary peaks supports texture-component redistribution.
Peak broadening or splitting supports increasing orientation heterogeneity. None
of these results alone identifies a specific active slip or twin system.

## 5. Module 11: spatial localization of intragranular gradients

### 5.1 Scientific question

This module determines where the orientation-gradient signal occurs:

- broadly within grain interiors;
- preferentially near low- or high-angle boundaries;
- in localized bands or connected regions;
- near unindexed areas or map edges, where the signal may be less reliable.

### 5.2 Registered KAM and GROD definitions

The existing condition-level KAM definition remains:

- first neighbour order;
- 5 degree neighbour-misorientation exclusion threshold;
- indexed Ti-Hex pixels;
- native, non-interpolated grid.

It is retained for compatibility with the existing summary. It is not the
primary metric for the boundary-distance profile because it may include
2–5 degree neighbour pairs and therefore be mathematically coupled to the LAGB
network being tested.

The primary spatial-profile metric is `intragranular_KAM`:

- first-order orthogonal four-neighbour topology;
- neighbours restricted to the same grain under the 2 degree reconstruction;
- cross-boundary, cross-phase, unindexed, and missing neighbours excluded;
- native, non-interpolated orientations.

Sensitivity cases use:

- first and second neighbour orders;
- 2 and 5 degree exclusion thresholds.

GROD uses the current grain-reference orientation definition and the registered
2 degree grain reconstruction with `minPixel = 5`. Because grain reconstruction
changes the grain reference and can subdivide a high-gradient grain, GROD
sensitivity uses 1, 2, and 5 degree reconstruction thresholds. Focused 1.5 and
2.5 degree cases are added if the 1–5 degree sensitivity changes the
interpretation.

For every state, report:

- valid-pixel count and fraction;
- KAM mean, median, P75, P90, and P95;
- GROD mean, median, P75, P90, and P95;
- area fractions with KAM greater than 0.5, 1.0, and 1.5 degrees.

### 5.3 Boundary-distance profiles

Construct separate native-grid masks for:

- 2–5 degree boundaries;
- 5–15 degree boundaries;
- boundaries at least 15 degrees.

Boundary geometry comes from the unsmoothed primary reconstruction. Pixel-centre
distance is calculated to the nearest geometric boundary segment, not to a
boundary-adjacent pixel centre. A spatial index or equivalent bounded search may
be used for performance, but the reported distance is point-to-segment distance.

For distances from 0 to 10 micrometres in 0.5 micrometre bins, report:

- eligible and valid pixel counts;
- KAM mean, median, P90, and high-KAM fraction;
- GROD mean, median, and P90;
- normalized value relative to the grain-interior reference band, where
  sufficient interior pixels exist.

Profiles are not interpreted when a bin has inadequate support. The exact count
is always written so that sparsely populated outer-distance bins cannot appear
equivalent to well-sampled bins.

Bin edge inclusivity is fixed in the parameter contract. A 1 micrometre
aggregation is reported as a sensitivity analysis because the primary 0.5
micrometre bin width equals the scan step. Each bin also reports its fraction of
the eligible mapped area and the fraction adjacent to unindexed or non-Ti-Hex
sites.

Because absolute boundary distance covaries with grain size and with increasing
boundary-network density, a secondary profile uses distance normalized by the
equivalent grain radius, and/or stratifies the profile by registered grain-size
classes. This secondary analysis is not allowed to replace the physical-distance
profile.

### 5.4 Localized-region analysis

For the three registered high-KAM thresholds, identify connected regions using
the native map topology. Report:

- region count;
- total area fraction;
- median and maximum region area;
- major-axis angle to AD;
- fraction touching unindexed regions or the map edge;
- fraction intersecting each boundary-distance class.

Region summaries distinguish extended localization from isolated high-value
pixels. The minimum-region-area rule is declared in the parameter file and
tested for sensitivity.

### 5.5 Masks and artefact control

- Unindexed pixels are never assigned KAM, GROD, or distance-profile values.
- Map-edge and unindexed-neighbour flags are retained.
- Primary KAM excludes neighbour pairs above the registered threshold.
- No gap filling is introduced for primary quantitative results.
- Smoothed or denoised maps may improve visual continuity but cannot replace
  original-data values.
- Changes in valid-pixel fraction are plotted or tabulated alongside metric
  changes.

### 5.6 Required outputs

```text
11_spatial_localization/
  spatial_localization_parameters.csv
  spatial_support_audit.csv
  spatial_gradient_summary.csv
  boundary_distance_profiles.csv
  high_gradient_regions.csv
  spatial_sensitivity_summary.csv
  spatial_gradient_trends.png
  boundary_distance_profiles.png
  spatial_gradient_maps.png
```

### 5.7 Interpretation boundary

Increasing KAM or GROD near a boundary class supports increasing lattice
curvature or orientation-gradient localization near that network. It does not
directly measure total dislocation density. A high-gradient region coincident
with unindexed areas is treated cautiously rather than automatically interpreted
as a deformation band. Pixel-level p-values are prohibited because spatially
correlated pixels are not independent experimental replicates.

## 6. Module 12: topology of low-angle boundary networks

### 6.1 Scientific question

Length fraction alone does not distinguish isolated low-angle segments from a
connected subboundary network. The primary connectivity network is the combined
2–15 degree LAGB population. This avoids introducing artificial breaks where the
local misorientation of one physical boundary crosses the 5 degree reporting
cutoff. The 2–5 and 5–15 degree subnetworks remain separate outputs for
describing the composition and apparent maturation of the LAGB population.

### 6.2 Graph construction

Use the unsmoothed MTEX grain-boundary geometry:

- geometric boundary vertices are graph nodes;
- retained boundary segments are weighted graph edges;
- edge weights are physical segment lengths;
- segment orientation is measured as an undirected acute angle to AD.

MTEX `ebsdId` endpoint pairs identify indexed pixels on the sides of a boundary
segment and are retained for provenance. They are not substituted for geometric
boundary vertices.

Coincident geometric vertices are merged under a tolerance substantially smaller
than the 0.5 micrometre scan step. Duplicate edges and zero-length segments are
audited. Metrics are computed from unsmoothed geometry; smoothed geometry is
display-only.

Native MTEX boundary segments are pixel-scale faces and do not constitute
physically meaningful branch-length observations. After graph construction,
degree-2 chains are compressed into branches between endpoints and junctions.
Closed degree-2 cycles are retained as explicitly labelled cycle branches.
Length-distribution statistics use these branches, not individual pixel-scale
faces.

Junctions generated by raster stair-steps or ideal four-way pixel-grid
intersections are tested using vertex-clustering radii of 0, 0.5, and 1
micrometre. Short-branch pruning at 0, 1, and 2 micrometres is a mandatory
sensitivity test. Grain-detection threshold and `minPixel` sensitivity are
registered because the 2 degree network is especially sensitive to
reconstruction.

### 6.3 Network metrics

For each boundary class and scan, report:

- segment count and total length;
- physical length density;
- branch-length median, P75, P90, and maximum;
- node density;
- node-degree distribution;
- endpoint density;
- junction density for degree at least 3;
- connected-component count and density;
- isolated-component length fraction;
- largest-component length and its fraction of total network length;
- largest-component x and y span;
- largest-component weighted graph diameter, defined as the maximum finite
  shortest-path distance and never called the graph's “longest path”;
- length-weighted segment-angle distribution relative to AD;
- length fractions within 15 and 30 degrees of AD and transverse to AD.

Every component reports whether it touches a map edge or an unindexed/non-Ti-Hex
hole. Endpoint density is reported both as observed and after excluding these
censored endpoints.

### 6.4 Map-spanning criterion

An AD-spanning network must contact both registered left and right edge bands.
A transverse-spanning network must contact both bottom and top edge bands. The
edge-band width is declared in pixels and subjected to sensitivity testing.

If a component does not meet the explicit criterion, the manuscript uses
“increased connectivity”, “larger connected component”, or equivalent wording,
not “percolation” or “through-going network”. Even when the criterion is met, it
is described as a two-dimensional field-of-view-spanning network, not a
three-dimensional percolating boundary network.

### 6.5 Required outputs

```text
12_lagb_topology/
  lagb_topology_parameters.csv
  lagb_topology_summary.csv
  lagb_components.csv
  lagb_branches.csv
  lagb_topology_sensitivity.csv
  lagb_topology_trends.png
  lagb_network_maps.png
```

### 6.6 Interpretation boundary

An increase in low-angle-boundary length and network connectivity is consistent
with increasing organization of orientation gradients into subboundary
structures. Conventional EBSD alone does not prove that every detected segment
is a mature dislocation cell wall, nor does it establish the underlying
dislocation character. Segment-angle results describe two-dimensional boundary
traces in the mapped section and are not interpreted as three-dimensional grain-
boundary-plane orientations.

## 7. Module 13: two-dimensional EBSD minimum-energy GND estimate

### 7.1 Scientific question and formal name

The module evaluates whether the observable in-plane lattice-curvature signal
changes with reduction and whether it localizes spatially. Every table, figure,
and manuscript reference uses the formal name:

> minimum-energy GND estimate constrained by two-dimensional EBSD observable
> lattice curvature

The result may be described as a lower-bound-type, model-dependent estimate, but
not as a strict mathematical lower bound, a complete Nye-tensor result, or total
dislocation density.

### 7.2 Curvature and Nye-tensor limitation

The native indexed Ti-Hex orientation field is assigned the registered `grainId`,
gridified, and used to calculate orientation curvature. In a two-dimensional
EBSD map, MTEX evaluates the x and y derivatives and leaves the third, surface-
normal derivative column unavailable. The curvature tensor and the corresponding
Nye/dislocation-density tensor are therefore incomplete; only six curvature
components enter the fit.

The missing derivative is not imputed as measured information. The fitted value
represents a minimum-energy solution compatible with the observable curvature
components and the assumed dislocation-system basis.

### 7.3 HCP titanium dislocation-system basis

The candidate basis includes, where supported by the installed MTEX version:

- basal `<a>` edge and screw;
- prismatic `<a>` edge and screw;
- pyramidal `<a>` edge and screw;
- first-order pyramidal `<c+a>` edge and screw;
- second-order pyramidal `<c+a>` edge and screw.

Three registered candidate sets are compared:

- `a_only`;
- `a_plus_ca1`;
- `a_plus_ca1_ca2`, used as the complete primary candidate set.

Burgers-vector magnitudes are derived from the registered Ti-Hex lattice
parameters in the CTF/MTEX crystal symmetry and converted explicitly to metres.
The current twelve CTF files consistently declare `a = 2.954` Angstrom and
`c = 4.729` Angstrom; the generator verifies this rather than silently assuming
it.

The exact crystallographic systems, multiplicities, line directions, Burgers
vectors, and relative line-energy assumptions are written to the parameter
table. Systems are constructed explicitly from the MTEX basal, prismatic-A,
pyramidal-A, pyramidal-CA, and pyramidal-2CA slip-system definitions, followed by
antipodal symmetrization and edge/screw dislocation construction. Duplicate
systems, multiplicities, and the resulting `<a>` and `<c+a>` Burgers-vector
lengths are asserted. No additional manual one-third scaling is applied after
the MTEX dislocation-system constructor.

The primary fit uses MTEX dislocation-system tensor construction and
`fitDislocationSystems` or its version-equivalent official method. No custom
inversion may silently replace the registered official formulation.

MTEX 6.1.1 does not provide a complete predefined HCP slip-system factory for
this purpose. The listed titanium systems must therefore be constructed
explicitly and tested; the workflow must not rely on
`dislocationSystem.hcp(cs)` to provide the candidate basis.

The runtime preflight verifies that MATLAB Optimization Toolbox `linprog` is
available and that the algorithm invoked by MTEX 6.1.1 works in the installed
MATLAB release. Fit failures remain `NaN` and contribute to an explicit success
fraction.

Two energy models are mandatory:

- equal weights, `u = 1`, as the transparent reference;
- a normalized isotropic line-energy approximation, with
  `u_edge = (b/b_a)^2` and
  `u_screw = (1 - nu) (b/b_a)^2`.

The CP-Ti Poisson ratio used in the second model requires a primary literature
source and a registered sensitivity range. Fitted system-family contributions
remain model-sensitivity outputs and are not interpreted as direct slip-system
activity.

### 7.4 Primary summaries

For every scan, report:

- total Ti-Hex pixels;
- pixels with valid observable curvature;
- fit-success pixels and fraction;
- valid fractions for each of the six observable curvature components and all
  exclusion reasons;
- fit reconstruction-residual median, P90, and maximum;
- median, P25, P75, P90, and P95 of
  `rho_l1_m2 = factor * sum(abs(rho), 2)`;
- median and P90 of the separately named normalized
  `energy_weighted_index_m2 = factor * sum(abs(rho) .* u, 2)`;
- mean and standard deviation of `log10(rho_l1_m2)` over finite positive values;
- distance-to-boundary profiles compatible with Module 11;
- raw-versus-denoised differences.

Zero, missing, and non-finite values are treated separately. Logarithms are not
computed by silently adding an arbitrary constant.

No single absolute "high GND" threshold is introduced without an independent
physical basis. Condition comparisons emphasize the median, interquartile
range, P90, common log-scale distributions, and profile shifts.

### 7.5 Boundary and data-quality handling

Before curvature is calculated, the indexed Ti-Hex data receive grain IDs from
the primary 2 degree reconstruction. Primary finite differences must not bridge:

- unindexed sites;
- phase boundaries;
- map gaps;
- grain boundaries excluded by the declared misorientation rule.

The analysis records boundary-exclusion width, valid stencil count, and
edge-pixel loss. A map with reduced valid coverage is not compared to another map
without reporting the coverage difference.

The primary result excludes at least one native step, 0.5 micrometres, from
boundaries, gaps, non-Ti-Hex sites, and map edges. Exclusion widths of 0, 0.5,
and 1.0 micrometres are sensitivity cases. Grain-reconstruction thresholds of
2, 2.5, and 5 degrees test whether automatic cross-grain derivative suppression
changes the result.

For `denoised_raw_common`, the mask additionally excludes phase-changed or
originally invalid coordinates and their one-pixel neighbourhood. Denoising-
induced coverage gains are never interpreted as physical GND gains or losses.

### 7.6 Mandatory sensitivity cases

The sensitivity table includes, at minimum:

- original versus denoised CTF;
- original unsmoothed data versus one declared within-grain smoothing case that
  preserves the original valid mask;
- alternative curvature neighbourhood or derivative stencil supported by MTEX;
- at least two boundary-exclusion widths;
- full registered HCP system set versus a documented reduced set;
- declared relative edge/screw line-energy assumptions.

Step size is fixed by the experiment and is not treated as an independent
experimental variable. A controlled 0.5 versus 1.0 micrometre within-grain
spatial-scale test is labelled a resolution sensitivity check, not new physical
evidence. It must not average or interpolate across 2 degree grain boundaries,
phase boundaries, or invalid sites.

The fit uses the conversion `factor` returned by `fitDislocationSystems`; the
workflow does not hard-code a multiplier. The output records scan length unit,
Ti-Hex lattice parameters, Burgers-vector lengths, and `factor`. For the current
micrometre/Angstrom inputs, a test verifies the expected order and the specific
returned conversion factor before a full run.

### 7.7 Required outputs and storage

```text
13_gnd_lower_bound/
  gnd_lower_bound_parameters.csv
  gnd_dislocation_systems.csv
  gnd_lower_bound_summary.csv
  gnd_lower_bound_distributions.csv
  gnd_boundary_distance_profiles.csv
  gnd_sensitivity.csv
  gnd_lower_bound_trends.png
  gnd_lower_bound_maps.png
  intermediate/
    per_state_compact_results.mat
```

The CSV distribution output stores histogram or quantile summaries rather than
millions of per-pixel rows. Compact MAT files preserve the arrays required to
reproduce maps and aggregate statistics. This avoids unnecessary duplication of
the already large result tree.

### 7.8 Interpretation boundary

The module can support statements that the observable in-plane curvature or its
minimum-energy GND estimate increases, decreases, or localizes. It cannot by
itself:

- quantify total mobile plus statistically stored dislocation density;
- prove a unique active slip system;
- prove that GND changes are the sole source of strengthening;
- establish a temporal sequence from one field per deformation condition;
- validate absolute density without an independent method.

TEM, X-ray diffraction line-profile analysis, or another independent defect
characterization method is required before making strong absolute-dislocation
density claims.

## 8. Cross-module synthesis

The four modules are interpreted jointly using the following decision structure.

### 8.1 Evidence compatible with early net c-axis reorientation

This interpretation requires a tracked texture peak or component centroid to
move crystallographically, with a corresponding change in c-axis direction.
A change in mean c-axis angle without stable peak tracking is insufficient.

### 8.2 Evidence compatible with texture-component redistribution

This interpretation is preferred when peak locations are nearly stationary but
their integrated fractions or MRD values exchange. It is distinct from a
continuous rotation mechanism.

### 8.3 Evidence compatible with increasing substructure organization

This interpretation is strengthened when several of the following co-occur:

- increasing KAM/GROD or high-gradient area fraction;
- stronger boundary-proximal localization;
- increasing 2–5 or 5–15 degree boundary length density;
- increasing largest-component fraction or junction density;
- increasing observable-curvature-constrained minimum-energy GND estimate.

No single threshold metric is treated as decisive.

### 8.4 Evidence compatible with heterogeneous or non-monotonic evolution

If texture, gradient, topology, or GND quantities change non-monotonically, the
paper reports the measured sequence and evaluates raw/denoised and parameter
sensitivity. A smooth two-stage mechanism is not imposed on inconsistent data.

## 9. Relationship to tensile data

Condition-level original-data metrics may be joined to the accepted tensile
summary. The join may examine descriptive correspondence with:

- `Rp0.2`;
- `Rm`;
- yield ratio;
- elongation;
- reduction of area;
- reliable uniform-strain and work-hardening quantities.

With six deformation conditions:

- Spearman correlations are descriptive;
- leave-one-condition-out ranges are mandatory;
- pixel or grain counts are not experimental replication;
- denoised scans are not additional observations;
- causal claims require mechanism-specific corroboration.

Texture peak motion, boundary-network connectivity, and minimum-energy GND metrics
will not be selected after inspecting tensile correlations solely to maximize
correlation strength.

## 10. Experimental supplementation guided by outcomes

The analysis determines which additional experiments have the highest value:

- The first priority for a high-level paper is at least three independent EBSD
  fields per deformation state, with a registered sampling design spanning the
  bar centre, mid-radius, and near-surface region where specimen geometry
  permits. The present single 300 by 300 micrometre field per state supports a
  descriptive mapped-field sequence but not material-wide variance estimates.
- If a texture peak migrates but the active deformation mode remains ambiguous,
  use deformation-path modelling, in-situ or interrupted diffraction, or
  slip-trace/TEM evidence.
- If KAM/GROD and LAGB topology increase but GND absolute values are unstable,
  use TEM dislocation/subboundary imaging and XRD line-profile analysis.
- If boundary-localized gradients dominate, acquire replicated EBSD fields at
  different radial and circumferential positions and consider higher-resolution
  EBSD.
- If peaks split or spatial bands appear, use higher-resolution maps,
  serial-section/3D methods where justified, and targeted TEM lift-outs.
- If mechanical trends are non-monotonic, add tensile repeats and interrupted
  deformation states before assigning a stage transition.

## 11. Validation strategy

### 11.1 Pure synthetic tests

Create deterministic tests for:

- crystallographic peak matching under symmetry;
- peak birth, disappearance, split, merge, and ambiguous cases;
- boundary-distance bins on a known raster geometry;
- high-KAM connected-region area and orientation;
- graph component, node degree, spanning, and weighted-diameter metrics;
- Burgers-vector and density unit conversion;
- GND fit output shape, sign handling, and finite-value masks.

### 11.2 One-pair formal smoke test

Before all twelve inputs are processed, run one original/denoised pair through
all four modules and verify:

- input hashes are unchanged;
- AD is horizontal;
- expected files and columns exist;
- figures are visually legible;
- raw and denoised rows are paired;
- validity and coverage counts are plausible;
- compact intermediates reproduce plotted summaries.

### 11.3 Full integration checks

The full run must verify:

- exact coverage of six samples times two variants;
- registered sample order and deformation values;
- parameter-table coverage for every run;
- summary fractions and counts remain in physical ranges;
- topology edge lengths conserve the retained segment length;
- peak fractions and track identities are internally consistent;
- no finite difference bridges an excluded GND stencil;
- output-contract inventory is complete;
- source CTF SHA-256 values are unchanged;
- MATLAB `checkcode`, module tests, integration tests, image inspection, and
  `git diff --check` pass.

## 12. Failure handling

- A missing or invalid output causes an explicit failed module status; it is not
  replaced by an empty-looking figure.
- A peak that cannot be matched receives an ambiguous or unmatched status.
- A sparse distance bin retains its sample count and is excluded from
  interpretation.
- A network without a spanning component is reported as non-spanning.
- A GND fit with insufficient valid curvature coverage is flagged and omitted
  from trend interpretation rather than silently compared.
- Raw/denoised disagreement is reported as a sensitivity result.

## 13. Parallel implementation boundaries

After this design and the subsequent implementation plan are approved, parallel
work may proceed with disjoint ownership:

1. texture-component detection, tracking, and tests;
2. spatial-localization and LAGB-topology analysis, with separately owned files;
3. minimum-energy GND feasibility, implementation, and tests;
4. main-agent ownership of shared contract, manifest, orchestration, integration,
   and final validation.

All parallel agents use `gpt-5.6-sol` with high reasoning effort. Shared files are
modified only by the main agent after module interfaces are frozen. Shared files
include the output contract and its tests, main runner and integration tests,
manifest writer, README/inventory logic, generic raw/denoised comparison, and
tensile integration.

## 14. Completion criterion

The extension is complete when:

- all four modules process the twelve registered scans;
- original-data conclusions and denoised-data sensitivity are separated;
- every output is registered in the contract and manifest;
- automated and visual validation passes;
- the four modules can distinguish peak motion, component redistribution,
  orientation-gradient localization, boundary-network organization, and a
  method-dependent minimum-energy GND estimate;
- the final interpretation states whether the existing data support a staged,
  overlapping, non-monotonic, or unresolved microstructural evolution;
- unsupported absolute dislocation-density and unique-mechanism claims are
  excluded.
