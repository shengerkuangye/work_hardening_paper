# Ti-Hex Grain-Boundary Misorientation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract spatially connected Ti-Hex/Ti-Hex boundary segments from six raw EBSD scans and compare their length-weighted misorientation distributions and low-angle-boundary metrics.

**Architecture:** Two pure MATLAB helpers perform connected-component length filtering and weighted angle summarization. One MTEX batch function reconstructs boundaries at three detection floors, combines the extracted inner and outer Ti-Hex/Ti-Hex segments, classifies every segment from its measured misorientation, validates length conservation, writes tables, and renders statistical comparison figures.

**Tech Stack:** MATLAB R2025a, MTEX 6.1.1, MATLAB tables and graphics, CTF EBSD input.

## Global Constraints

- Use only the six original CTF files; do not use or overwrite `_denoised.ctf` files.
- Import with `convertEuler2SpatialReferenceFrame`.
- Reconstruct with detection floors 0.5, 1, and 2 degrees and a 15-degree high-angle threshold.
- Use `minPixel = 5`, do not smooth boundaries, and retain only Ti-Hex/Ti-Hex segments.
- Filter inner-boundary connected components shorter than 1.0 micrometre from primary statistics while retaining unfiltered audit totals.
- Classify LAGB and HAGB from each segment's measured angle, not from the MTEX object name.
- Weight distributions and fractions by `segLength`.
- Write derived files only under `results/mtex_grain_boundary_misorientation/`.
- Do not produce plan-view EBSD maps.

---

## File Structure

- Create `tools/mtex/component_length_mask.m`: pure connected-component length filter.
- Create `tools/mtex/summarize_weighted_boundary_angles.m`: pure weighted histogram and descriptive statistics.
- Create `tools/mtex/generate_grain_boundary_misorientation_distribution.m`: MTEX extraction, sensitivity analysis, tables, validation, and plots.
- Create `tools/mtex/test_grain_boundary_misorientation_distribution.m`: synthetic unit tests and full-data integration assertions.
- Create five derived files under `results/mtex_grain_boundary_misorientation/` as specified in the design.

### Task 1: Pure boundary filtering and weighted distribution helpers

**Files:**
- Create: `tools/mtex/component_length_mask.m`
- Create: `tools/mtex/summarize_weighted_boundary_angles.m`
- Create: `tools/mtex/test_grain_boundary_misorientation_distribution.m`

**Interfaces:**
- Consumes component identifiers, segment lengths, and minimum component length; produces a logical segment-retention mask and per-component total length.
- Consumes angle and segment-length vectors plus fixed bin edges; produces a weighted statistics struct and distribution table.
- Produces `[keepMask,componentLength] = component_length_mask(componentId,segLength,minLength)`.
- Produces `[stats,distribution] = summarize_weighted_boundary_angles(thetaDeg,segLengthUm,binEdgesDeg)`.

- [ ] **Step 1: Write failing synthetic tests**

Create the test function with these assertions before either helper exists:

```matlab
function test_grain_boundary_misorientation_distribution(scanRoot, outputDir)
arguments
  scanRoot (1,1) string = ""
  outputDir (1,1) string = ""
end

componentId = [1;1;2;3;3];
segLengthUm = [0.4;0.7;0.9;0.25;0.25];
[keepMask, componentLength] = component_length_mask( ...
  componentId, segLengthUm, 1.0);
assert(isequal(keepMask, logical([1;1;0;0;0])));
assert(isequal(componentLength, [1.1;0.9;0.5]));
[emptyKeep, emptyLength] = component_length_mask([], [], 1.0);
assert(isempty(emptyKeep) && isempty(emptyLength));

thetaDeg = [1.2;2.2;15.2];
weights = [1;3;2];
edges = 1:0.5:94;
[stats, distribution] = summarize_weighted_boundary_angles( ...
  thetaDeg, weights, edges);
assert(stats.segment_count == 3);
assert(abs(stats.total_boundary_length_um - 6) < 1e-12);
assert(abs(stats.weighted_mean_angle_deg - 38.2/6) < 1e-12);
assert(abs(stats.weighted_median_angle_deg - 2.2) < 1e-12);
assert(abs(sum(distribution.length_fraction) - 1) < 1e-12);
assert(abs(sum(distribution.probability_density_per_degree .* ...
  distribution.bin_width_deg) - 1) < 1e-12);

if scanRoot ~= ""
  assert(outputDir ~= "", "outputDir is required for integration testing.");
  [summary, sensitivity, fullDistribution] = ...
    generate_grain_boundary_misorientation_distribution(scanRoot, outputDir);
  assert(height(summary) == 6);
  assert(height(sensitivity) == 18);
  assert(height(fullDistribution) == 1116);
  expectedSamples = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"];
  assert(isequal(summary.sample, expectedSamples));
  assert(all(summary.total_retained_boundary_length_um > 0));
  assert(all(abs(summary.lagb_length_um + summary.hagb_length_um - ...
    summary.total_retained_boundary_length_um) < 1e-8));
  intervalSum = summary.length_1_2_um + summary.length_2_5_um + ...
    summary.length_5_10_um + summary.length_10_15_um + ...
    summary.length_15_94_um;
  assert(all(abs(intervalSum - summary.total_retained_boundary_length_um) ...
    < 1e-8));
  assert(isequal(unique(sensitivity.detection_floor_deg), [0.5;1;2]));
  for sample = expectedSamples'
    rows = fullDistribution.sample == sample;
    assert(nnz(rows) == 186);
    assert(abs(sum(fullDistribution.length_fraction(rows)) - 1) < 1e-10);
  end
  expectedFiles = [
    "grain_boundary_misorientation_distribution.csv";
    "grain_boundary_misorientation_summary.csv";
    "grain_boundary_detection_sensitivity.csv";
    "grain_boundary_misorientation_distribution.png";
    "grain_boundary_misorientation_metrics.png"
  ];
  assert(all(isfile(fullfile(outputDir, expectedFiles))));
end

fprintf("GRAIN_BOUNDARY_MISORIENTATION_TESTS_OK\n");
end
```

- [ ] **Step 2: Run the test and verify RED**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath('C:/Users/22069/Documents/GitHub/work_hardening_paper/tools/mtex'); test_grain_boundary_misorientation_distribution"
```

Expected: nonzero exit because `component_length_mask` is undefined.

- [ ] **Step 3: Implement the component-length filter**

Create `component_length_mask.m`:

```matlab
function [keepMask, componentLength] = component_length_mask( ...
  componentId, segLength, minLength)
componentId = double(componentId(:));
segLength = double(segLength(:));
assert(numel(componentId) == numel(segLength));
if isempty(componentId)
  keepMask = false(0,1);
  componentLength = zeros(0,1);
  return
end
assert(all(componentId >= 1 & componentId == fix(componentId)));
assert(all(isfinite(segLength) & segLength > 0));
assert(isscalar(minLength) && isfinite(minLength) && minLength >= 0);
componentLength = accumarray(componentId, segLength, [], @sum);
keepMask = componentLength(componentId) >= minLength;
keepMask = logical(keepMask(:));
end
```

- [ ] **Step 4: Implement weighted statistics**

Create `summarize_weighted_boundary_angles.m`. Validate equal nonempty finite vectors, positive lengths, angles inside the requested edges, and exact `1:0.5:94` edges. Compute:

```matlab
weightedMean = sum(thetaDeg .* segLengthUm) / sum(segLengthUm);
[sortedTheta, order] = sort(thetaDeg);
sortedWeight = segLengthUm(order);
weightedMedian = sortedTheta(find(cumsum(sortedWeight) >= ...
  0.5 * sum(sortedWeight), 1));
[boundaryLengthUm,~,binIndex] = histcounts(thetaDeg, binEdgesDeg);
boundaryLengthUm = accumarray(binIndex(binIndex > 0), ...
  segLengthUm(binIndex > 0), [numel(binEdgesDeg)-1,1], @sum, 0);
segmentCount = accumarray(binIndex(binIndex > 0), 1, ...
  [numel(binEdgesDeg)-1,1], @sum, 0);
```

Return statistics fields `segment_count`, `total_boundary_length_um`,
`weighted_mean_angle_deg`, `weighted_median_angle_deg`, `min_angle_deg`,
`max_angle_deg`, and `mode_bin_center_deg`. Return distribution fields
`bin_lower_deg`, `bin_upper_deg`, `bin_width_deg`, `bin_center_deg`,
`segment_count`, `boundary_length_um`, `length_fraction`,
`probability_density_per_degree`, and `cumulative_length_fraction`.

- [ ] **Step 5: Run the synthetic test and verify GREEN**

Run the Step 2 command. Expected: exit 0 and
`GRAIN_BOUNDARY_MISORIENTATION_TESTS_OK`.

- [ ] **Step 6: Commit Task 1**

```powershell
git add -- tools/mtex/component_length_mask.m tools/mtex/summarize_weighted_boundary_angles.m tools/mtex/test_grain_boundary_misorientation_distribution.m
git commit -m "feat: summarize length-weighted grain boundaries"
```

### Task 2: MTEX boundary extraction and sensitivity tables

**Files:**
- Create: `tools/mtex/generate_grain_boundary_misorientation_distribution.m`
- Use: `tools/mtex/component_length_mask.m`
- Use: `tools/mtex/summarize_weighted_boundary_angles.m`
- Test: `tools/mtex/test_grain_boundary_misorientation_distribution.m`

**Interfaces:**
- Consumes scalar string `scanRoot` and `outputDir`.
- Produces `[summary,sensitivity,distribution] = generate_grain_boundary_misorientation_distribution(scanRoot,outputDir)` and five derived files.

- [ ] **Step 1: Run the full-data test and verify integration RED**

Run the Task 1 test with the two paths from its interface. Expected: nonzero
exit because `generate_grain_boundary_misorientation_distribution` is
undefined.

- [ ] **Step 2: Implement six-scan, three-floor extraction**

Use the sample mapping in the design and the constants:

```matlab
detectionFloorsDeg = [0.5;1;2];
classificationAngleDeg = 15;
minPixel = 5;
minInnerComponentLengthUm = 1.0;
distributionEdgesDeg = 1:0.5:94;
```

For every scan, load once, calculate indexed Ti-Hex area as
`length(ebsd('Ti-Hex')) * polyarea(ebsd.unitCell.x,ebsd.unitCell.y)`, then
repeat reconstruction for each floor:

```matlab
[grains, ebsd.grainId] = calcGrains(ebsd('indexed'), ...
  'threshold', [floorDeg classificationAngleDeg] * degree, ...
  'minPixel', minPixel);
inner = grains.innerBoundary('Ti-Hex','Ti-Hex');
outer = grains.boundary('Ti-Hex','Ti-Hex');
innerTheta = double(angle(inner.misorientation) / degree);
outerTheta = double(angle(outer.misorientation) / degree);
innerTheta = innerTheta(:);
outerTheta = outerTheta(:);
innerLength = double(inner.segLength(:));
outerLength = double(outer.segLength(:));
[keepInner,~] = component_length_mask( ...
  inner.componentId, innerLength, minInnerComponentLengthUm);
thetaDeg = [innerTheta(keepInner); outerTheta];
segLengthUm = [innerLength(keepInner); outerLength];
valid = isfinite(thetaDeg) & isfinite(segLengthUm) & segLengthUm > 0 & ...
  thetaDeg >= floorDeg - 1e-8 & thetaDeg <= 94;
thetaDeg = thetaDeg(valid);
segLengthUm = segLengthUm(valid);
```

Record raw and retained inner lengths, outer length, removed length, grain
count, indexed Ti-Hex pixels, indexed Ti-Hex area, and the final measured
angle/length vectors. Never infer LAGB/HAGB membership from `inner` or
`outer`; use `thetaDeg < 15` versus `thetaDeg >= 15`.

- [ ] **Step 3: Build primary and sensitivity rows**

For every floor calculate total retained length, LAGB/HAGB lengths and
fractions, LAGB/HAGB density per indexed Ti-Hex area, length-weighted mean
and median, and audit values. Append one row to the 18-row sensitivity
table.

For the 1-degree floor, additionally compute exact interval masks:

```matlab
interval1 = thetaDeg >= 1 & thetaDeg < 2;
interval2 = thetaDeg >= 2 & thetaDeg < 5;
interval3 = thetaDeg >= 5 & thetaDeg < 10;
interval4 = thetaDeg >= 10 & thetaDeg < 15;
interval5 = thetaDeg >= 15 & thetaDeg <= 94;
```

Store each interval's length and fraction in the six-row primary summary.
Call `summarize_weighted_boundary_angles` with `1:0.5:94`, prepend sample
metadata, and append its 186 rows to the distribution table.

- [ ] **Step 4: Validate tables and write CSV files**

Assert row counts 6, 18, and 1116. For each sample and floor assert positive
indexed area and retained length, length conservation, fractions in `[0,1]`,
and finite angles. For the primary summary assert the five interval lengths
sum to total retained length. For every distribution block assert its length
fraction and density integral are one.

Write:

```matlab
writetable(summary, fullfile(outputDir, ...
  'grain_boundary_misorientation_summary.csv'));
writetable(sensitivity, fullfile(outputDir, ...
  'grain_boundary_detection_sensitivity.csv'));
writetable(distribution, fullfile(outputDir, ...
  'grain_boundary_misorientation_distribution.csv'));
```

- [ ] **Step 5: Run integration test and verify tables GREEN**

Run the full-data command. Expected: all six scans run at all three floors,
the three CSV files are written, and the test advances to the expected
missing-figure assertion.

### Task 3: Statistical comparison figures and final verification

**Files:**
- Modify: `tools/mtex/generate_grain_boundary_misorientation_distribution.m`
- Create: `results/mtex_grain_boundary_misorientation/grain_boundary_misorientation_distribution.png`
- Create: `results/mtex_grain_boundary_misorientation/grain_boundary_misorientation_metrics.png`
- Test: `tools/mtex/test_grain_boundary_misorientation_distribution.m`

**Interfaces:**
- Consumes the primary summary and distribution tables from Task 2.
- Produces two 300-dpi white-background PNG files without axes toolbars.

- [ ] **Step 1: Implement the distribution figure**

Create an invisible 1200-by-900 figure with two tiled axes. Plot the six
length-weighted probability-density curves with consistent sample colors.
Use `[1,94]` for the complete panel and `[1,15]` for the low-angle detail.
Label the x axis `Ti-Hex grain-boundary misorientation (degree)` and the y
axis `Length-weighted probability density (degree^{-1})`.

- [ ] **Step 2: Implement the metrics figure**

Create an invisible 1200-by-800 two-panel figure. Plot LAGB length fraction
in percent against cold reduction in the upper panel and LAGB length density
in micrometre per square micrometre in the lower panel. Label every point
with its sample name and keep the cold-reduction ordering from 7d to 5d.

- [ ] **Step 3: Export and canonicalize both figures**

Remove axes toolbars, call `exportgraphics(...,'Resolution',300,
'BackgroundColor','white')`, then round-trip each file through
`imread`/`imwrite` to remove volatile PNG metadata.

- [ ] **Step 4: Run full integration test and verify GREEN**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath('C:/Users/22069/Documents/GitHub/work_hardening_paper/tools/mtex'); test_grain_boundary_misorientation_distribution('C:/Users/22069/Documents/GitHub/work_hardening_paper/data/ebsd_kpl_250221_7_df/scans','C:/Users/22069/Documents/GitHub/work_hardening_paper/results/mtex_grain_boundary_misorientation')"
```

Expected: exit 0 and `GRAIN_BOUNDARY_MISORIENTATION_TESTS_OK`.

- [ ] **Step 5: Run static and structural checks**

Run MATLAB `checkcode(...,'-id')` on all four new scripts and assert no
issues. Run `git diff --check`. Read all three CSVs back and repeat row-count,
length-conservation, and probability-conservation assertions.

- [ ] **Step 6: Inspect both PNG files at original resolution**

Verify that all six curves are visible, the 1-15 degree detail is legible,
metric labels are not clipped, and no plan-view EBSD map was produced.

- [ ] **Step 7: Commit the verified analysis**

```powershell
git add -- tools/mtex/generate_grain_boundary_misorientation_distribution.m tools/mtex/test_grain_boundary_misorientation_distribution.m results/mtex_grain_boundary_misorientation
git commit -m "feat: compare Ti-Hex grain-boundary misorientations"
```
