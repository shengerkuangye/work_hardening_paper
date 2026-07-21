# Comprehensive EBSD analysis from original and denoised CTF data

## 1. Objective

Build a reproducible MTEX workflow that uses only the six original CTF scans and their six denoised counterparts. The workflow will replace vendor-generated EBSD charts with consistently reconstructed maps, quantitative tables, raw-versus-denoised sensitivity analyses, and manuscript-ready figures.

The horizontal map direction is the bar axial direction (AD), corresponding to the MTEX specimen x direction. The vertical map direction is the transverse/radial in-plane direction, corresponding to the MTEX specimen y direction. The scan normal is the MTEX specimen z direction.

The primary scientific questions are:

1. How do grain morphology, crystallographic texture, boundary populations, and intragranular orientation gradients evolve with cold reduction?
2. Which trends are robust to denoising, reconstruction threshold, and angular cutoff choices?
3. Do the data support net c-axis reorientation, texture-component redistribution, early dislocation-substructure formation, later substructure maturation, or spatial heterogeneity?
4. Which EBSD quantities can be related cautiously to the tensile response without treating six condition-level observations as independent mechanistic proof?

## 2. Input data and sample order

Use the following six scan folders and their paired `*_denoised.ctf` files under `data/ebsd_kpl_250221_7_df/scans/`:

| Sample | Folder | Original CTF | Denoised CTF | Cold reduction (%) |
| --- | --- | --- | --- | ---: |
| 7d | `d7` | `ebsd_sample_7_map_15.ctf` | `ebsd_sample_7_map_15_denoised.ctf` | 0.00 |
| 6.48d | `d6_48` | `ebsd_sample_648_map_13.ctf` | `ebsd_sample_648_map_13_denoised.ctf` | 14.31 |
| 6.02d | `d6_02` | `ebsd_sample_602_map_11.ctf` | `ebsd_sample_602_map_11_denoised.ctf` | 26.04 |
| 5.6d | `d5_6` | `ebsd_sample_56_map_9.ctf` | `ebsd_sample_56_map_9_denoised.ctf` | 36.00 |
| 5.25d | `d5_25` | `ebsd_sample_525_map_7.ctf` | `ebsd_sample_525_map_7_denoised.ctf` | 43.75 |
| 5d | `d5` | `ebsd_sample_5_map_3.ctf` | `ebsd_sample_5_map_3_denoised.ctf` | 48.98 |

All scans are expected to retain the native 600 by 600 grid at 0.5 micrometre steps. Raw CTF files are read-only.

## 3. Evidentiary role of original and denoised data

- Original CTF data are the primary source for quantitative conclusions.
- Denoised CTF data are processed with exactly the same workflow for paired sensitivity analysis and, when appropriate, clearer visual maps.
- A denoised map may be used in a manuscript figure only when its caption identifies the processing and the corresponding quantitative trend is consistent with the original data.
- If denoising materially changes a metric, the original value remains primary and the difference is reported rather than hidden.
- Original and denoised scans are not treated as independent experimental replicates.

## 4. Common import, geometry, and reconstruction rules

- Import CTF orientations using `convertEuler2SpatialReferenceFrame`, consistent with the validated project workflows.
- Preserve the complete native grid, including unindexed sites; do not interpolate or bridge missing sites for primary calculations.
- Restrict crystallographic calculations to indexed Ti-Hex pixels unless a specific audit reports the cubic-Ti or unindexed fraction.
- Retain persistent EBSD IDs when mapping boundary endpoints to native coordinates.
- Use `unitCell` grain reconstruction so unindexed regions do not create nonlocal adjacency.
- Apply no boundary smoothing in primary quantitative calculations. Smoothed boundary maps may be produced only as explicitly labelled visual variants.
- Use a 2 degree detection floor, 15 degree HAGB threshold, and `minPixel = 5` as the primary reconstruction. Sensitivity cases will test relevant floors and thresholds.
- Boundary statistics are weighted by physical segment length; count-weighted results are secondary sensitivity outputs.
- All map figures use the same sample orientation: AD horizontal, positive x to the right.

## 5. Analysis modules

### 5.1 Data-quality and raw-denoised audit

For every scan pair, report:

- total, indexed, Ti-Hex, cubic-Ti, and unindexed pixel counts and fractions;
- grid dimensions, step sizes, mapped area, duplicate coordinates, and missing coordinate combinations;
- MAD, BC, BS, band-count, and error-code distributions where available;
- pixels whose phase, index status, or orientation changed during denoising;
- symmetry-reduced orientation-change distributions for paired indexed Ti-Hex pixels;
- spatial maps of filled, removed, phase-changed, and materially reoriented pixels;
- SHA-256 hashes of all twelve inputs before and after analysis.

### 5.2 Standard map atlas

Generate consistently scaled maps for original and denoised data:

- IPF-X/AD, IPF-Y, and IPF-Z;
- band-contrast, MAD, phase, unindexed-pixel, and grain-ID maps;
- IPF and BC maps with boundary overlays for 2-5, 5-15, and at least 15 degrees;
- candidate alpha-Ti twin-boundary overlays when crystallographic matching is valid;
- individual maps and six-condition composite panels.

### 5.3 Grain morphology

Report distributions and condition summaries for:

- grain area, perimeter, and equivalent-circle diameter;
- area-weighted equivalent diameter;
- ellipse major/minor axes and aspect ratio;
- maximum/minimum Feret diameters where supported robustly;
- circularity, elongation, and shape factor;
- long-axis angle relative to AD;
- grain-shape orientation tensor or equivalent axial-alignment metric;
- the fraction of grains exceeding predefined, visible elongation thresholds.

Equivalent diameter will not be used alone to infer refinement. Morphology metrics must distinguish grain subdivision from directional elongation.

### 5.4 Boundary and misorientation analysis

Produce:

- length-weighted boundary-misorientation distributions;
- count-weighted distributions as sensitivity results;
- boundary lengths, fractions, and indexed-area-normalized densities for 2-5, 5-10, 10-15, and at least 15 degrees;
- detection-floor sensitivity at 0.5, 1, and 2 degrees;
- direct 5 versus 15 degree LAGB upper-cutoff comparisons on the same boundary population;
- misorientation-axis distributions where statistically and crystallographically meaningful;
- candidate twin-boundary length and area fractions with stated angular tolerances;
- raw-versus-denoised boundary-network differences.

Existing validated native-grid topology checks and conservation tests will be reused rather than duplicated inconsistently.

### 5.5 Intragranular deformation metrics

Calculate and compare:

- KAM maps and distributions with stated neighbour order and exclusion cutoff;
- neighbour-order and cutoff sensitivity;
- grain-reference orientation deviation angle and axis (GROD);
- grain orientation spread (GOS);
- within-grain orientation-range or equivalent orientation-dispersion metrics;
- grain-level relations among KAM/GROD/GOS, size, aspect ratio, and boundary proximity;
- Nye-tensor/GND estimates if MTEX input quality and available curvature components permit a defensible calculation.

Any GND result is reported as a step-size- and method-dependent lower-bound proxy for geometrically necessary dislocations, not total dislocation density. Raw-versus-denoised sensitivity is mandatory for KAM, GROD, GOS, and GND-related results.

### 5.6 Texture and c-axis analysis

Recalculate texture directly from the CTF orientations:

- `{0001}`, `{10-10}`, and `{11-20}` pole figures;
- AD, in-plane transverse/radial, and scan-normal inverse pole figures;
- ODF sections using one documented kernel and halfwidth, with sensitivity where necessary;
- maximum MRD, texture J-index, M-index, and an orientation-concentration/entropy metric where well defined;
- major fibre/component volume fractions with explicit angular tolerances;
- c-axis angles to AD, in-plane transverse/radial direction, and scan normal;
- full c-axis spherical distributions and azimuthal distributions around AD;
- pole-figure peak positions, widths, and multimodality diagnostics;
- pairwise ODF-distance or texture-change metrics between adjacent deformation states;
- raw-versus-denoised texture comparisons.

The analysis must distinguish net c-axis tilt from azimuthal rotation, texture spreading, peak splitting, and exchange among texture components.

### 5.7 Axial-tensile orientation factors

For subsequent uniaxial tension along AD, calculate maps and distributions for:

- basal `<a>` slip;
- prismatic `<a>` slip;
- pyramidal `<a>` slip;
- pyramidal `<c+a>` slip;
- relevant alpha-Ti twin systems when supported;
- maximum Schmid factor by family;
- an orientation/Taylor-factor estimate using a documented constitutive assumption.

These outputs describe orientation favourability for the later axial tensile test. They must not be presented as direct identification of active rotary-swaging mechanisms without a validated swaging stress path or deformation model.

### 5.8 EBSD-tensile integration

Create a condition-level table joining the primary original-data EBSD metrics to the accepted tensile summary:

- `Rp0.2`, `Rm`, yield ratio, elongation, reduction of area, and reliable uniform-strain metrics;
- grain morphology and axial-alignment metrics;
- LAGB/HAGB fractions and length densities;
- KAM/GROD/GOS and any defensible GND proxy;
- texture indices, component fractions, and c-axis metrics;
- axial-tensile Schmid/Taylor metrics.

With only six deformation conditions, correlations are descriptive and hypothesis-generating. Pixel counts are never treated as experimental replication, and causal language is prohibited.

## 6. Output structure

Write new derived artifacts under `results/mtex_ebsd_comprehensive/`:

```text
00_audit/
01_standard_maps/
02_grain_morphology/
03_boundary_analysis/
04_intragranular_deformation/
05_texture/
06_schmid_taylor/
07_raw_denoised_comparison/
08_tensile_correlation/
README.md
analysis_manifest.csv
```

Each quantitative module produces auditable CSV files, individual figures, six-condition composite figures, and paired raw-versus-denoised comparison outputs where relevant. Figure export uses consistent sample order, colours, axis limits, units, labels, and publication resolution.

## 7. Implementation organization

The workflow will be divided into shared MTEX utilities and bounded module generators. Common sample metadata, import, coordinate definitions, native-grid audits, reconstruction parameters, figure styling, and output-manifest logic will be centralized so modules cannot silently diverge.

Parallel implementation may be used for independent modules after shared interfaces are fixed:

1. audit and raw-denoised comparison;
2. grain morphology plus standard maps;
3. texture plus c-axis analysis;
4. boundary, intragranular, Schmid/Taylor, and tensile integration in staged follow-up work.

All parallel agents use `gpt-5.6-sol` with `high` reasoning effort. Agents share the repository, so parallel tasks must own disjoint files or coordinate through the main agent before editing shared utilities.

## 8. Validation

Automated and integration checks must verify:

- all twelve inputs exist, retain expected hashes, and remain unchanged;
- every scan contains the native 600 by 600 coordinate grid at 0.5 micrometre steps;
- boundary endpoint pairs are native face neighbours;
- boundary interval lengths conserve eligible totals and fractions sum to one;
- grain and pixel statistics remain within valid ranges;
- orientation angles and Schmid factors remain within physical bounds;
- paired raw-denoised rows exist for every sample and metric;
- map axes and AD annotations are consistent with AD = specimen x;
- synthetic orientation tests validate c-axis, boundary-classification, morphology, KAM/GROD, and Schmid-factor calculations;
- MATLAB `checkcode`, module tests, full twelve-scan integration, expected-output inventory, image inspection, and `git diff --check` pass before delivery.

## 9. Interpretation boundaries

- LAGB and KAM are indirect evidence of orientation gradients and dislocation organization, not direct total-dislocation-density measurements.
- GND estimates do not include statistically stored dislocations and depend on angular precision, step size, and smoothing.
- A single scan per condition cannot establish material-wide spatial variance; reported condition trends remain scan-representative until replicated radial/field measurements are available.
- EBSD grain reconstruction does not by itself demonstrate physical recrystallization, recovery, or a unique slip/twin mechanism.
- Equivalent grain diameter is not sufficient evidence of refinement when grains become elongated.
- Texture and Schmid-factor calculations establish orientation tendencies, not active-system proof.
- Mechanism claims must be separated into measured results, supported inferences, and hypotheses requiring TEM, diffraction, replicated EBSD, or modelling.

## 10. Completion criterion

The comprehensive workflow is complete when all twelve CTF files run through the same reproducible pipeline, every agreed module produces validated CSV and figure outputs, raw-versus-denoised sensitivity is explicit, the original inputs remain unchanged, and the resulting README documents which metrics are suitable for manuscript conclusions and which remain provisional.
