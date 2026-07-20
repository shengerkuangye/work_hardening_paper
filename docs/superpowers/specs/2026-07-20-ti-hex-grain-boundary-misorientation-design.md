# Ti-Hex Boundary-Segment Misorientation Distribution Design

## Objective

Compare the length-weighted misorientation distributions of thresholded
Ti-Hex/Ti-Hex boundary segments in the six original EBSD scans. The analysis
tests how the measured low-angle segment fraction and length density vary with
cold reduction and with the angular detection floor.

The estimator is an MTEX thresholded Ti-Hex/Ti-Hex boundary-segment network.
It is not a connected-wall length measurement, and no plan-view boundary map
is produced.

## Input data and native geometry

Use only the original CTF files under
`data/ebsd_kpl_250221_7_df/scans/`:

| Sample | Folder | CTF file | Cold reduction (%) |
| --- | --- | --- | ---: |
| 7d | d7 | ebsd_sample_7_map_15.ctf | 0.00 |
| 6.48d | d6_48 | ebsd_sample_648_map_13.ctf | 14.31 |
| 6.02d | d6_02 | ebsd_sample_602_map_11.ctf | 26.04 |
| 5.6d | d5_6 | ebsd_sample_56_map_9.ctf | 36.00 |
| 5.25d | d5_25 | ebsd_sample_525_map_7.ctf | 43.75 |
| 5d | d5 | ebsd_sample_5_map_3.ctf | 48.98 |

Import with `convertEuler2SpatialReferenceFrame`. Each imported object must
retain the complete native 600 by 600 grid, including unindexed sites. The
required audit values are 360000 mapped sites, 0.5 micrometre steps in both
directions, `scanUnit = "um"`, 0.25 square micrometre per pixel, and 90000
square micrometres of mapped area.

Do not use the denoised CTF files. Do not subset to indexed sites, gridify,
interpolate, or otherwise infill the native grid.

## Grain reconstruction and topology

For each detection floor `delta`, reconstruct from the full EBSD object:

```matlab
[grains, ebsd.grainId] = calcGrains(ebsd, 'unitCell', ...
  'threshold', [delta 15] * degree, 'minPixel', 5);
```

The `unitCell` option preserves native face adjacency in the presence of
unindexed regions without bridging missing sites.
For every Ti-Hex/Ti-Hex boundary segment, map `grainBoundary.ebsdId` through
the persistent IDs in `ebsd.id`; do not interpret those IDs as table row
indices. The mapped endpoints must be horizontal or vertical nearest
neighbours on the 0.5 micrometre native grid. Reject the reconstruction if any
Ti-Hex/Ti-Hex source face is nonlocal.

No boundary smoothing or connected-component length post-filter is applied.

## Source-by-angle operational definition

The boundary source and the measured symmetry-reduced misorientation jointly
define eligibility. At each detection floor `delta`, partition all finite,
positive-length Ti-Hex/Ti-Hex source segments into four mutually exclusive
cells:

| MTEX source | Measured angle | Role |
| --- | --- | --- |
| `grains.innerBoundary` | `delta <= theta < 15 degrees` | eligible LAGB |
| `grains.innerBoundary` | `theta >= 15 degrees` | audit-only inner high-angle |
| `grains.boundary` | `delta <= theta < 15 degrees` | audit-only outer low-angle |
| `grains.boundary` | `theta >= 15 degrees` | eligible HAGB |

The analysis population is exactly the union of eligible inner LAGB and
eligible outer HAGB segments. Inner high-angle segments must not enter the HAGB
population, and outer low-angle segments must not enter the LAGB population.
Record counts and lengths for all four cells and verify both source-level and
eligible/excluded conservation.

## Detection-floor sensitivity and statistics

- Detection floors: 0.5, 1, and 2 degrees.
- Primary detection floor: 1 degree.
- LAGB/HAGB classification angle: 15 degrees.
- Grain reconstruction minimum: 5 pixels.
- Boundary smoothing: none.
- Segment weighting: MTEX `segLength`, never segment count.

For all three floors, report eligible LAGB/HAGB segment counts, lengths,
fractions, length densities normalized by indexed Ti-Hex area, and
length-weighted mean and median angle. The metrics figure must plot the 0.5-,
1-, and 2-degree series for both LAGB fraction and LAGB density so that floor
sensitivity is visible.

For the primary 1-degree result, use unsmoothed 0.5-degree bins over 1 to 94
degrees and also report the exact intervals 1-2, 2-5, 5-10, 10-15, and 15-94
degrees. The distribution figure must identify 1 degree as the primary
detection floor.

Because the low-angle network changes materially with the detection floor,
trend interpretation must compare all three sensitivity series. The 1-degree
result is the designated primary estimate, not a claim of detection-floor
invariance.

## Audit fields

The sensitivity and primary summary tables record:

- scan unit, native dimensions and steps, pixel area, mapped count and area;
- unindexed count, indexed Ti-Hex count and area, and indexed Ti-Hex fraction
  of mapped area;
- Ti-Hex grain count only;
- counts and lengths for both MTEX sources and all four source-by-angle cells;
- exhaustive source-network, excluded-cross-class, and eligible totals;
- persistent-ID endpoint-pair count, nonlocal-pair count, and maximum endpoint
  distance;
- eligible LAGB/HAGB metrics and length-weighted angle statistics.

The six-row primary summary must exactly reproduce the 1-degree rows of the
18-row sensitivity table before adding the ten interval columns. The primary
distribution must conserve the corresponding summary segment count and
boundary length for every sample.

## Outputs

Write exactly five derived artifacts under
`results/mtex_grain_boundary_misorientation/`:

- `grain_boundary_misorientation_summary.csv`;
- `grain_boundary_detection_sensitivity.csv`;
- `grain_boundary_misorientation_distribution.csv`;
- `grain_boundary_misorientation_distribution.png`;
- `grain_boundary_misorientation_metrics.png`.

The outputs contain statistical comparisons only. No plan-view EBSD map is
part of this analysis.

## Validation

The implementation must verify that:

- every scan is the full native 600 by 600 Cartesian grid with unique
  coordinates at 0.5 micrometre steps, all 360000 site combinations, and a
  positive unindexed count;
- every source segment has Ti-Hex on both sides, a positive length, a finite
  symmetry-reduced angle, and native face-adjacent endpoints;
- all four source-by-angle cells are exhaustive and exclusive;
- eligible LAGB plus eligible HAGB equals the complete eligible population;
- interval lengths and fractions conserve the primary eligible total;
- every length-weighted distribution integrates to one and reproduces its
  primary summary row;
- all six samples and all three sensitivity floors are present;
- the output directory contains the exact five freshly generated artifacts;
- raw CTF hashes are unchanged after generation.
