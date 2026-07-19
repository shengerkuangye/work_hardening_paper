# Ti-Hex Grain-Boundary Misorientation Distribution Design

## Objective

Extract spatially connected Ti-Hex/Ti-Hex grain-boundary segments from the
six original EBSD CTF scans, then compare their length-weighted
misorientation distributions. The analysis is intended to test whether the
amount of low-angle grain boundary increases with cold deformation.

The output is a statistical grain-boundary analysis. It will not include a
plan-view grain-boundary map.

## Input Data

Use the six original CTF files under
`data/ebsd_kpl_250221_7_df/scans/`:

| Sample | Folder | CTF file | Cold reduction (%) |
| --- | --- | --- | ---: |
| 7d | d7 | ebsd_sample_7_map_15.ctf | 0.00 |
| 6.48d | d6_48 | ebsd_sample_648_map_13.ctf | 14.31 |
| 6.02d | d6_02 | ebsd_sample_602_map_11.ctf | 26.04 |
| 5.6d | d5_6 | ebsd_sample_56_map_9.ctf | 36.00 |
| 5.25d | d5_25 | ebsd_sample_525_map_7.ctf | 43.75 |
| 5d | d5 | ebsd_sample_5_map_3.ctf | 48.98 |

Do not use or overwrite the `_denoised.ctf` files. Import using
`convertEuler2SpatialReferenceFrame`.

The original CTF files do not contain grain identifiers or independent
grain-boundary labels. Their vendor-generated Mackenzie distributions omit
angles below approximately 5 degrees and therefore cannot answer the
low-angle-boundary question.

## Boundary Definition

The analysis distinguishes the boundary detection floor from the later
low-angle/high-angle classification.

- Primary detection floor: 1 degree.
- Sensitivity detection floors: 0.5 and 2 degrees.
- Low-angle/high-angle classification: 15 degrees.
- Low-angle grain boundary (LAGB): detected Ti-Hex/Ti-Hex boundary with
  misorientation from the detection floor up to, but not including,
  15 degrees.
- High-angle grain boundary (HAGB): detected Ti-Hex/Ti-Hex boundary with
  misorientation of at least 15 degrees.

For each detection floor `delta`, reconstruct grains using MTEX dual
thresholds `[delta, 15] * degree`. MTEX stores the boundaries from `delta`
to 15 degrees as `grains.innerBoundary` and boundaries of at least 15
degrees as `grains.boundary`. `innerBoundary` is an extracted subgrain
boundary; it is not treated as an ordinary intragranular pixel pair.

Use only boundary segments with Ti-Hex on both sides and finite
misorientation. Exclude scan edges, phase boundaries, and boundaries next
to unindexed pixels. Use `minPixel = 5` during reconstruction to suppress
isolated indexed speckles. Do not smooth the reconstructed boundaries
before measuring their length.

Suppress a boundary component from the primary statistics when its total
connected length is less than 1.0 micrometre. Record the removed segment
count and length so that this continuity filter remains auditable. The
unfiltered totals will also be retained in the sensitivity table.

## Statistics

Weight every accepted boundary segment by its MTEX `segLength`. Do not
weight by segment count or by the number of neighboring pixels.

For the primary 1-degree detection floor, use common bins from 1 to 94
degrees with a 0.5-degree width. For each sample, report:

- total detected Ti-Hex/Ti-Hex boundary length;
- LAGB length and HAGB length;
- LAGB fraction of total detected boundary length;
- HAGB fraction of total detected boundary length;
- boundary-length fractions in 1-2, 2-5, 5-10, 10-15, and 15-94 degree
  intervals;
- LAGB and HAGB length density normalized by indexed Ti-Hex scan area;
- length-weighted mean and length-weighted median misorientation angle.

Repeat the principal length and fraction metrics for the 0.5-, 1-, and
2-degree detection floors. A deformation trend will be treated as robust
only if its direction is retained across these sensitivity cases.

## Outputs

Write new derived files under `results/mtex_grain_boundary_misorientation/`:

- `grain_boundary_misorientation_distribution.csv`: the primary
  length-weighted 1-94 degree distributions for all six samples;
- `grain_boundary_misorientation_summary.csv`: primary per-sample lengths,
  fractions, densities, and descriptive statistics;
- `grain_boundary_detection_sensitivity.csv`: metrics at 0.5-, 1-, and
  2-degree detection floors, including unfiltered and continuity-filtered
  boundary length;
- `grain_boundary_misorientation_distribution.png`: two panels showing the
  complete 1-94 degree distribution and the 1-15 degree detail;
- `grain_boundary_misorientation_metrics.png`: LAGB fraction and LAGB
  length density versus cold reduction.

The plots compare samples only; no EBSD plan-view map will be generated.

## Validation

The implementation must verify that:

- every retained segment has Ti-Hex on both sides, positive length, and a
  finite symmetry-reduced misorientation within the HCP range;
- the integrated length-weighted probability of each distribution is one;
- LAGB length plus HAGB length equals the total retained boundary length;
- interval lengths sum to the total retained boundary length;
- reported densities use positive indexed Ti-Hex area;
- all six samples and all three detection floors are present;
- raw input files remain unchanged.

The final interpretation will distinguish measured changes in boundary
length/fraction from possible effects of indexing loss, angular resolution,
and the selected detection floor.
