# LAGB 5° versus 15° Threshold Comparison Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate reproducible MTEX statistics and figures comparing 2°–5° and 2°–15° Ti-Hex low-angle boundary definitions for all six EBSD scans.

**Architecture:** Add one pure MATLAB metric helper and one independent batch generator that reuses the validated native-grid boundary partition and weighted-histogram helpers. Keep the detection floor fixed at 2°, reconstruct one common `[2°,15°]` boundary population from the full EBSD grid, classify that same population at 5° and 15°, and write all new artifacts to a separate results directory.

**Tech Stack:** MATLAB R2025a, MTEX 6.1.1, native MATLAB tables/graphics, PowerShell verification.

## Global Constraints

- Read only the six original CTF scans; never use `_denoised.ctf` inputs.
- Keep the detection floor fixed at 2°, reconstruct once with `[2°,15°]`, and
  compare post-reconstruction classification angles 5° and 15° on the same
  eligible boundary population.
- Use full-grid `calcGrains(...,'unitCell',...,'minPixel',5)` without smoothing.
- Analyze indexed Ti-Hex/Ti-Hex boundary segments and weight statistics by segment length.
- Do not overwrite `results/mtex_grain_boundary_misorientation/`.
- Describe results as an MTEX-thresholded boundary-segment network, not manually traced physical walls.

---

### Task 1: Pure threshold metric helper

**Files:**
- Create: `tools/mtex/calculate_boundary_threshold_metrics.m`
- Create: `tools/mtex/test_grain_boundary_threshold_comparison.m`

**Interfaces:**
- Consumes: angle vector in degrees, positive boundary-segment lengths in µm, detection floor, classification angle, and indexed Ti-Hex area in µm².
- Produces: `stats = calculate_boundary_threshold_metrics(thetaDeg, segLengthUm, detectionFloorDeg, classificationAngleDeg, indexedAreaUm2)`.

- [ ] **Step 1: Write the failing synthetic test**

Add a test that calls the nonexistent helper for angles
`[2.1;4.9;5;10;14.9;15;30]` and lengths `[1;2;3;4;5;6;7]`. Assert that the
5° case has LAGB length 3, the 15° case has LAGB length 15, each case conserves
total length 28, and the 5° LAGB length does not exceed the 15° value.

- [ ] **Step 2: Run the test and verify RED**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath('tools/mtex'); test_grain_boundary_threshold_comparison"
```

Expected: failure because `calculate_boundary_threshold_metrics` is undefined.

- [ ] **Step 3: Implement the minimal helper**

The helper must validate finite angles, positive lengths,
`classificationAngleDeg > detectionFloorDeg`, and positive indexed area. Use
`thetaDeg < classificationAngleDeg` for LAGB and its complement for HAGB, then
return segment counts, lengths, fractions, total density, LAGB density, HAGB
density, and length-weighted mean/median angle.

- [ ] **Step 4: Run the test and verify GREEN**

Run the Step 2 command. Expected: `GRAIN_BOUNDARY_THRESHOLD_COMPARISON_TESTS_OK`.

- [ ] **Step 5: Commit the tested helper**

```powershell
git add -- tools/mtex/calculate_boundary_threshold_metrics.m tools/mtex/test_grain_boundary_threshold_comparison.m
git commit -m "feat: calculate boundary threshold metrics"
```

### Task 2: Six-scan MTEX comparison generator

**Files:**
- Create: `tools/mtex/generate_grain_boundary_threshold_comparison.m`
- Modify: `tools/mtex/test_grain_boundary_threshold_comparison.m`
- Create: `results/mtex_grain_boundary_threshold_comparison/grain_boundary_threshold_comparison_summary.csv`
- Create: `results/mtex_grain_boundary_threshold_comparison/grain_boundary_threshold_comparison_distribution.csv`

**Interfaces:**
- Consumes: `scanRoot` and `outputDir` string scalars.
- Produces: `[summary, distribution] = generate_grain_boundary_threshold_comparison(scanRoot, outputDir)` and the five named result artifacts.

- [ ] **Step 1: Extend the test before implementation**

Require exact output names, 12 summary rows, 2,232 distribution rows, sample
order `7d` through `5d`, classification angles `[5;15]`, fixed detection floor
2°, zero nonlocal native-grid pairs, exact fraction conservation, and fresh
artifacts only. Add a static assertion prohibiting indexed-only reconstruction,
`gridify`, interpolation, and smoothing.

- [ ] **Step 2: Run the integration test and verify RED**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath('tools/mtex'); test_grain_boundary_threshold_comparison('data/ebsd_kpl_250221_7_df/scans','results/mtex_grain_boundary_threshold_comparison')"
```

Expected: failure because the generator is undefined.

- [ ] **Step 3: Implement the batch generator**

For each scan, load the original CTF, audit the 600×600 native grid, and call
once:

```matlab
[grains, ebsdFull.grainId] = calcGrains(ebsdFull, 'unitCell', ...
  'threshold', [2 15] * degree, 'minPixel', 5);
```

Select Ti-Hex/Ti-Hex `innerBoundary` and `boundary`, pass their angles and
lengths through `partition_ti_hex_boundary_segments` using the conventional
15° reconstruction classification, and verify persistent EBSD endpoint
adjacency with `audit_native_grid_pairs`. Apply the Task 1 helper to this same
angle/length population first at 5° and then at 15°, and summarize
`1:0.5:94` bins with
`summarize_weighted_boundary_angles`. Write the two CSV files after validating
all rows and conservation identities.

- [ ] **Step 4: Run the integration test and verify GREEN**

Run the Step 2 command. Expected: all six scans complete and the test prints
`GRAIN_BOUNDARY_THRESHOLD_COMPARISON_TESTS_OK`.

- [ ] **Step 5: Commit generator and tabular outputs**

```powershell
git add -- tools/mtex/generate_grain_boundary_threshold_comparison.m tools/mtex/test_grain_boundary_threshold_comparison.m results/mtex_grain_boundary_threshold_comparison/*.csv
git commit -m "feat: compare 5 and 15 degree boundary thresholds"
```

### Task 3: Matched figures and final verification

**Files:**
- Modify: `tools/mtex/generate_grain_boundary_threshold_comparison.m`
- Modify: `tools/mtex/test_grain_boundary_threshold_comparison.m`
- Create: `results/mtex_grain_boundary_threshold_comparison/grain_boundary_threshold_5deg.png`
- Create: `results/mtex_grain_boundary_threshold_comparison/grain_boundary_threshold_15deg.png`
- Create: `results/mtex_grain_boundary_threshold_comparison/grain_boundary_threshold_comparison.png`

**Interfaces:**
- Consumes: validated summary and distribution tables from Task 2.
- Produces: two matched four-panel threshold figures and one two-panel direct metric comparison.

- [ ] **Step 1: Add failing figure assertions**

Require all three nonempty PNG files to be newer than the generation start and
require the generator source to contain titles identifying the fixed 2° floor
and the relevant 5° or 15° upper threshold.

- [ ] **Step 2: Run the test and verify RED**

Run the Task 2 integration command. Expected: failure because the three PNGs do
not yet exist.

- [ ] **Step 3: Implement plots**

For each upper threshold, create a 2×2 figure containing full misorientation
distribution, low-angle detail, LAGB length fraction, and LAGB length density.
Use identical sample colors and comparable axes. Create a separate 1×2 figure
that overlays the 5° and 15° length fraction and length density series using
both color and marker differences. Export at 300 dpi.

- [ ] **Step 4: Run full automated verification**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath('tools/mtex'); test_grain_boundary_threshold_comparison('data/ebsd_kpl_250221_7_df/scans','results/mtex_grain_boundary_threshold_comparison'); files={'tools/mtex/calculate_boundary_threshold_metrics.m','tools/mtex/generate_grain_boundary_threshold_comparison.m','tools/mtex/test_grain_boundary_threshold_comparison.m'}; for k=1:numel(files), assert(isempty(checkcode(files{k},'-id'))); end; disp('FINAL_THRESHOLD_COMPARISON_OK');"
git diff --check
```

Expected: both success markers, MATLAB exit code 0, and no Git whitespace errors.

- [ ] **Step 5: Inspect all three PNGs at original resolution**

Verify legible legends, unclipped labels, identical sample order, explicit 2°
detection floor, and clear 5°/15° classification labels. Revise and repeat Step
4 if any visual defect is present.

- [ ] **Step 6: Commit final figures**

```powershell
git add -- tools/mtex/generate_grain_boundary_threshold_comparison.m tools/mtex/test_grain_boundary_threshold_comparison.m results/mtex_grain_boundary_threshold_comparison/*.png
git commit -m "feat: plot LAGB threshold comparison"
```
