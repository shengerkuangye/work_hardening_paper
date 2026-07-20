# Ti-Hex Native-Grid Boundary-Segment Implementation Plan

> **For agentic workers:** use test-driven development for every behavior
> change and retain RED/GREEN evidence from both synthetic and six-scan tests.

**Goal:** Generate auditable, length-weighted Ti-Hex/Ti-Hex boundary-segment
distributions and detection-floor sensitivity metrics from the six original
EBSD scans without altering native topology.

**Architecture:** Two pure helpers partition measured boundary segments and
audit native endpoint adjacency. A third pure helper calculates weighted
histograms. One MTEX batch generator loads each full scan once, reconstructs
grains at three floors with `unitCell`, validates all conservation identities,
writes three exact-schema tables, and renders two statistical figures.

**Tech stack:** MATLAB R2025a, MTEX 6.1.1, CTF EBSD input, MATLAB tables and
graphics.

## Global constraints

- Use the six original CTF files only; never use `_denoised.ctf` inputs.
- Import with `convertEuler2SpatialReferenceFrame`.
- Preserve the full 600 by 600 EBSD object, including unindexed sites.
- Require 360000 mapped sites, 0.5 micrometre x/y steps, and scan unit `um`.
- Reconstruct exactly with `calcGrains(ebsd, 'unitCell', 'threshold',
  [floorDeg 15] * degree, 'minPixel', 5)`.
- Never use an indexed-only subset, `gridify`, interpolation, or smoothing.
- Map `grainBoundary.ebsdId` through persistent `ebsd.id` values; do not treat
  the IDs as row indices.
- Reject any nonlocal Ti-Hex/Ti-Hex source face.
- Apply no connected-component length filter.
- Use measured symmetry-reduced angle for the four source-by-angle cells.
- Include only inner LAGB plus outer HAGB segments in analysis distributions
  and metrics.
- Weight all distributions, fractions, and densities by `segLength`.
- Write only the five requested derived artifacts; generate no plan-view map.

## Final file structure

- `tools/mtex/audit_native_grid_pairs.m`
- `tools/mtex/partition_ti_hex_boundary_segments.m`
- `tools/mtex/summarize_weighted_boundary_angles.m`
- `tools/mtex/generate_grain_boundary_misorientation_distribution.m`
- `tools/mtex/test_grain_boundary_misorientation_distribution.m`
- five derived artifacts under `results/mtex_grain_boundary_misorientation/`

## Task 1: Specify the corrected behavior with failing tests

**Files:**

- Modify `tools/mtex/test_grain_boundary_misorientation_distribution.m`.

- [x] Add a shuffled persistent-ID example that succeeds only after IDs are
  mapped to coordinates.
- [x] Add nonlocal and unmapped endpoint rejection cases.
- [x] Add a synthetic four-cell partition with inner-low, inner-high,
  outer-low, and outer-high segments.
- [x] Assert exact count and length conservation for both sources, eligible
  cells, and excluded cross-class cells.
- [x] Assert exact schemas for sensitivity, summary, and distribution tables.
- [x] Assert the summary equals the 1-degree sensitivity rows exactly.
- [x] Assert each distribution reproduces its primary summary count and
  length, integrates to one, and carries matching metadata.
- [x] Assert scan geometry, unit, area, unindexed count, indexed Ti-Hex count
  and area, Ti-Hex grain count, and zero nonlocal source faces.
- [x] Assert all six samples have 0.5-, 1-, and 2-degree rows.
- [x] Clear the output folder before integration and require the exact five
  fresh artifacts afterward.

**RED command:**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath('C:/Users/22069/Documents/GitHub/work_hardening_paper/tools/mtex'); test_grain_boundary_misorientation_distribution('C:/Users/22069/Documents/GitHub/work_hardening_paper/data/ebsd_kpl_250221_7_df/scans','C:/Users/22069/Documents/GitHub/work_hardening_paper/results/mtex_grain_boundary_misorientation')"
```

Expected initial failure: the previous generator calls the removed continuity
filter after loading the first scan.

## Task 2: Implement persistent-ID native-grid auditing

**Files:**

- Create `tools/mtex/audit_native_grid_pairs.m`.

- [x] Accept an N by 2 endpoint-ID matrix, the full EBSD persistent-ID vector,
  and the corresponding native coordinates.
- [x] Resolve endpoints with an explicit ID lookup.
- [x] Require unique positive mapped IDs and reject missing endpoint IDs.
- [x] Reject duplicated coordinates, omitted Cartesian sites, and irregular
  native-axis spacing.
- [x] Require each pair to differ by exactly one native x or y step and zero
  steps in the other direction.
- [x] Return endpoint count, nonlocal count, and maximum distance.

**Focused GREEN command:** run the no-argument test function and require
`GRAIN_BOUNDARY_MISORIENTATION_TESTS_OK`.

## Task 3: Implement four-cell source-by-angle partitioning

**Files:**

- Create `tools/mtex/partition_ti_hex_boundary_segments.m`.

- [x] Validate finite symmetry-reduced angles and positive segment lengths.
- [x] Partition inner source segments into eligible LAGB and audit-only
  high-angle cells.
- [x] Partition outer source segments into audit-only low-angle and eligible
  HAGB cells.
- [x] Return only inner LAGB plus outer HAGB for downstream statistics.
- [x] Return counts and lengths for all four cells and assert exhaustive source
  and eligible/excluded conservation.

## Task 4: Replace the reconstruction and output path

**Files:**

- Modify `tools/mtex/generate_grain_boundary_misorientation_distribution.m`.

- [x] Load each original CTF once as a complete EBSD object.
- [x] Audit unit, 600 by 600 dimensions, 0.5 micrometre steps, pixel area,
  mapped area/count, unindexed count, and indexed Ti-Hex count/area/fraction.
- [x] Reconstruct all three floors from the full object with `unitCell`, dual
  `[floorDeg 15]` thresholds, and `minPixel = 5`.
- [x] Store Ti-Hex grain count only.
- [x] Extract Ti-Hex/Ti-Hex inner and outer sources and use their measured
  symmetry-reduced angles.
- [x] Audit native adjacency over the complete source network before applying
  source-by-angle eligibility.
- [x] Calculate length-weighted eligible LAGB/HAGB statistics and densities.
- [x] At the primary 1-degree floor, calculate 1-2, 2-5, 5-10, 10-15, and
  15-94 degree intervals plus unsmoothed 0.5-degree histogram bins.
- [x] Validate the exact 6-row, 18-row, and 1116-row table relationships before
  writing CSV files.

## Task 5: Render the two statistical figures

**Files:**

- Modify `tools/mtex/generate_grain_boundary_misorientation_distribution.m`.

- [x] Plot the primary 1-degree length-weighted distribution over 1-94 degrees
  and its 1-15 degree detail.
- [x] State the primary 1-degree detection floor and absence of smoothing in
  the distribution figure.
- [x] Plot all three floor series in both LAGB fraction and LAGB density panels.
- [x] Label the primary 1-degree sample points and identify that series in the
  legend.
- [x] Export both figures at 300 dpi on white backgrounds, remove axes
  toolbars, and canonicalize PNG metadata.
- [x] Produce no plan-view map.

## Task 6: Final verification and reporting

- [x] Run the complete six-scan integration command from Task 1 and require
  exit code 0.
- [x] Run MATLAB `checkcode(...,'-id')` on all five MATLAB files and require
  zero issues.
- [x] Read back all three CSVs and independently recheck exact row counts,
  source-by-angle conservation, eligible conservation, summary/sensitivity
  equality, distribution/summary equality, and all three floor series.
- [x] Search the implementation for forbidden indexed-only, gridification,
  interpolation, smoothing, and denoised-input paths.
- [x] Run `git diff --check`.
- [x] Compare SHA-256 hashes of all six raw CTF files before and after the run.
- [x] Inspect both PNGs at original resolution for complete curves, legible
  labels, explicit floor identification, and three-floor metric series.
- [x] Record numerical changes, RED/GREEN evidence, commands, image QA,
  self-review, and concerns in `.superpowers/sdd/final-fix-report.md`.
- [x] Commit with subject `fix: preserve native-grid grain-boundary topology`.
