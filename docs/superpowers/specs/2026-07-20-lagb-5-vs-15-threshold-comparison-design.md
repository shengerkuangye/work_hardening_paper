# Ti-Hex LAGB 5° versus 15° threshold comparison

## Objective

Compare how the interpretation of the same six Ti-Hex EBSD scans changes when
the upper limit of a low-angle boundary is set to 5° or 15°. The comparison is
intended to distinguish the very-low-angle 2°–5° network from the broader
2°–15° low-angle network associated with deformation substructure.

## Fixed analysis choices

- Input scans: the six original CTF files already used by the validated MTEX
  grain-boundary workflow (`7d`, `6.48d`, `6.02d`, `5.6d`, `5.25d`, `5d`).
- Phase: indexed Ti-Hex/Ti-Hex pairs only.
- Detection/noise floor: 2° for both comparisons.
- Classification case 1: LAGB is `2° <= theta < 5°`; HAGB is `theta >= 5°`.
- Classification case 2: LAGB is `2° <= theta < 15°`; HAGB is `theta >= 15°`.
- Grain reconstruction: full native EBSD grid, `unitCell`, `minPixel = 5`, no
  boundary smoothing, and the corresponding two MTEX thresholds.
- Weighting: boundary-segment length, not segment count.
- Raw CTF files remain read-only.

## Outputs

Write new derived files under
`results/mtex_grain_boundary_threshold_comparison/` without replacing the
existing grain-boundary results:

1. `grain_boundary_threshold_comparison_summary.csv`: one row per sample and
   classification threshold, containing LAGB/HAGB lengths, fractions, and
   length densities.
2. `grain_boundary_threshold_comparison_distribution.csv`: length-weighted
   0.5° histogram data for both threshold cases.
3. `grain_boundary_threshold_5deg.png`: six-sample distribution and metrics for
   the 5° upper limit.
4. `grain_boundary_threshold_15deg.png`: the same layout for the 15° upper
   limit.
5. `grain_boundary_threshold_comparison.png`: aligned small multiples directly
   comparing LAGB length fraction and length density for the two definitions.

All figures must state that the detection floor is 2° and must use identical
sample order, axes, colors, and units where a direct comparison is intended.

## Interpretation limits

- The 2°–5° result represents the very-low-angle part of the boundary network;
  it must not be described as the complete conventional LAGB population.
- The 2°–15° result is the conventional LAGB definition used for the primary
  paper comparison.
- An increase in either network is indirect EBSD evidence of orientation
  gradients, dislocation-wall development, or subgrain structure. It is not a
  direct measurement of total dislocation density or GND density.
- The reported lengths describe the MTEX-thresholded Ti-Hex/Ti-Hex
  boundary-segment network, not manually traced physical walls.

## Verification

- A synthetic unit test must fail before the new comparison implementation
  exists and then pass after implementation.
- For every sample and threshold, LAGB plus HAGB length must equal total
  eligible boundary length and their fractions must sum to one.
- The 5° LAGB network must be a subset of the 15° LAGB angle interval.
- All twelve sample-threshold rows and all expected histogram bins must be
  present.
- Full six-scan MATLAB integration, `checkcode`, PNG visual inspection, and
  `git diff --check` must pass before delivery.
