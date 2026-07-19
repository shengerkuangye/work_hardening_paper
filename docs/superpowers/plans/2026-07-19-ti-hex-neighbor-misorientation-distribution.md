# Ti-Hex Neighbor Misorientation Distribution Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute and compare the unthresholded four-neighbor Ti-Hex misorientation distributions for six EBSD scans so that a later grain-boundary detection threshold can be selected from the observed low-angle behavior.

**Architecture:** Two pure MATLAB helpers generate unique four-neighbor pairs and summarize angle vectors into fixed bins. One MTEX batch function gridifies each CTF, computes symmetry-reduced Ti-Hex neighbor misorientations, assembles per-sample and pooled CSV tables, and exports a two-panel comparison figure.

**Tech Stack:** MATLAB R2025a, MTEX 6.1.1, MATLAB tables/graphics, CTF EBSD input.

## Global Constraints

- Do not modify or overwrite raw CTF data.
- Use the six original CTF files and exclude `_denoised.ctf` files.
- Import with `convertEuler2SpatialReferenceFrame`.
- Include only pairs whose two pixels are indexed as `Ti-Hex`.
- Use unique first-order horizontal and vertical neighbors; exclude diagonals.
- Apply Ti-Hex crystal symmetry when computing the minimum misorientation.
- Use common 0°–94° bins with fixed 0.1° width.
- Do not smooth, interpolate, reconstruct grains, or apply a grain-boundary threshold.
- Write derived files under `results/mtex_neighbor_misorientation/`.

---

## File Structure

- Create `tools/mtex/four_neighbor_pairs.m`: pure logical-grid neighbor indexing.
- Create `tools/mtex/summarize_misorientation_angles.m`: pure histogram and descriptive statistics.
- Create `tools/mtex/generate_neighbor_misorientation_distribution.m`: six-scan MTEX analysis, pooling, CSV output, and plotting.
- Create `tools/mtex/test_neighbor_misorientation_distribution.m`: synthetic and full-data assertions.
- Create `results/mtex_neighbor_misorientation/neighbor_misorientation_distribution.csv`.
- Create `results/mtex_neighbor_misorientation/neighbor_misorientation_summary.csv`.
- Create `results/mtex_neighbor_misorientation/neighbor_misorientation_distribution.png`.

### Task 1: Pure neighbor indexing and distribution statistics

**Files:**
- Create: `tools/mtex/four_neighbor_pairs.m`
- Create: `tools/mtex/summarize_misorientation_angles.m`
- Create: `tools/mtex/test_neighbor_misorientation_distribution.m`

**Interfaces:**
- Consumes: a 2-D logical phase mask; returns two equal-length vectors of linear indices.
- Consumes: a finite angle vector and increasing bin edges; returns a scalar statistics struct and a distribution table.
- Produces: `[index1,index2] = four_neighbor_pairs(phaseMask)` and `[stats,distribution] = summarize_misorientation_angles(thetaDeg,binEdgesDeg)`.

- [ ] **Step 1: Write the failing synthetic tests**

Create `tools/mtex/test_neighbor_misorientation_distribution.m`:

```matlab
function test_neighbor_misorientation_distribution(scanRoot, outputDir)
arguments
  scanRoot (1,1) string = ""
  outputDir (1,1) string = ""
end

phaseMask = logical([1 1 0; 1 1 1]);
[index1, index2] = four_neighbor_pairs(phaseMask);
pairs = sort([index1, index2], 2);
expectedPairs = sort([1 3; 2 4; 4 6; 1 2; 3 4], 2);
assert(size(pairs,1) == 5);
assert(size(unique(pairs, "rows"),1) == 5);
assert(isequal(sortrows(pairs), sortrows(expectedPairs)));

binEdgesDeg = 0:0.1:94;
thetaDeg = [0.04; 0.14; 1.04; 93.84];
[stats, distribution] = summarize_misorientation_angles( ...
  thetaDeg, binEdgesDeg);
assert(stats.pair_count == 4);
assert(abs(stats.mean_angle_deg - mean(thetaDeg)) < 1e-12);
assert(abs(stats.mode_bin_center_deg - 0.05) < 1e-12);
assert(sum(distribution.pair_count) == 4);
assert(abs(sum(distribution.probability) - 1) < 1e-12);
assert(abs(distribution.cumulative_probability(end) - 1) < 1e-12);
assert(abs(sum(distribution.probability_density_per_degree .* ...
  distribution.bin_width_deg) - 1) < 1e-12);

if scanRoot ~= ""
  assert(outputDir ~= "", "outputDir is required for integration testing.");
  [summary, fullDistribution] = ...
    generate_neighbor_misorientation_distribution(scanRoot, outputDir);
  assert(height(summary) == 7);
  assert(height(fullDistribution) == 6580);
  assert(isequal(summary.sample, ...
    ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d";"pooled"]));
  assert(all(summary.pair_count(1:6) > 0 & ...
    summary.pair_count(1:6) <= 718800));
  assert(summary.pair_count(7) == sum(summary.pair_count(1:6)));
  for sample = summary.sample'
    rows = fullDistribution.sample == sample;
    assert(nnz(rows) == 940);
    assert(abs(sum(fullDistribution.probability(rows)) - 1) < 1e-10);
    assert(abs(fullDistribution.cumulative_probability(find(rows,1,"last")) ...
      - 1) < 1e-10);
  end
  assert(isfile(fullfile(outputDir, ...
    "neighbor_misorientation_distribution.csv")));
  assert(isfile(fullfile(outputDir, ...
    "neighbor_misorientation_summary.csv")));
  assert(isfile(fullfile(outputDir, ...
    "neighbor_misorientation_distribution.png")));
end

fprintf("NEIGHBOR_MISORIENTATION_TESTS_OK\n");
end
```

- [ ] **Step 2: Run the test and verify RED**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath('C:/Users/22069/Documents/GitHub/work_hardening_paper/tools/mtex'); test_neighbor_misorientation_distribution"
```

Expected: MATLAB exits nonzero because `four_neighbor_pairs` is undefined.

- [ ] **Step 3: Implement unique four-neighbor indexing**

Create `tools/mtex/four_neighbor_pairs.m`:

```matlab
function [index1, index2] = four_neighbor_pairs(phaseMask)
%FOUR_NEIGHBOR_PAIRS Return unique horizontal and vertical valid pairs.
assert(islogical(phaseMask) && ismatrix(phaseMask) && ...
  all(size(phaseMask) >= [2,2]), ...
  "phaseMask must be a logical matrix of at least 2-by-2.");
gridIndex = reshape(1:numel(phaseMask), size(phaseMask));

horizontalMask = phaseMask(:,1:end-1) & phaseMask(:,2:end);
horizontal1 = gridIndex(:,1:end-1);
horizontal2 = gridIndex(:,2:end);
verticalMask = phaseMask(1:end-1,:) & phaseMask(2:end,:);
vertical1 = gridIndex(1:end-1,:);
vertical2 = gridIndex(2:end,:);

index1 = [horizontal1(horizontalMask); vertical1(verticalMask)];
index2 = [horizontal2(horizontalMask); vertical2(verticalMask)];
index1 = index1(:);
index2 = index2(:);
end
```

- [ ] **Step 4: Implement fixed-bin statistics**

Create `tools/mtex/summarize_misorientation_angles.m` using population standard deviation, `prctile`, `histcounts`, and:

```matlab
thetaDeg = double(thetaDeg(:));
binEdgesDeg = double(binEdgesDeg(:)');
assert(~isempty(thetaDeg) && all(isfinite(thetaDeg)), ...
  "thetaDeg must contain finite values.");
assert(all(thetaDeg >= binEdgesDeg(1) & thetaDeg <= binEdgesDeg(end)), ...
  "Angles fall outside the requested bin range.");
assert(binEdgesDeg(1) == 0 && binEdgesDeg(end) == 94 && ...
  all(abs(diff(binEdgesDeg) - 0.1) < 1e-12), ...
  "binEdgesDeg must span 0 to 94 degrees in 0.1-degree bins.");

pairCount = histcounts(thetaDeg, binEdgesDeg)';
binLowerDeg = binEdgesDeg(1:end-1)';
binUpperDeg = binEdgesDeg(2:end)';
binWidthDeg = binUpperDeg - binLowerDeg;
binCenterDeg = (binLowerDeg + binUpperDeg) / 2;
probability = pairCount / numel(thetaDeg);
probabilityDensityPerDegree = probability ./ binWidthDeg;
cumulativeProbability = cumsum(probability);
[~, modeIndex] = max(pairCount);
quartiles = prctile(thetaDeg, [10,25,50,75,90]);

stats.pair_count = numel(thetaDeg);
stats.mean_angle_deg = mean(thetaDeg);
stats.std_angle_deg = std(thetaDeg,1);
stats.q10_angle_deg = quartiles(1);
stats.q25_angle_deg = quartiles(2);
stats.median_angle_deg = quartiles(3);
stats.q75_angle_deg = quartiles(4);
stats.q90_angle_deg = quartiles(5);
stats.mode_bin_center_deg = binCenterDeg(modeIndex);
stats.min_angle_deg = min(thetaDeg);
stats.max_angle_deg = max(thetaDeg);

distribution = table(binLowerDeg, binUpperDeg, binWidthDeg, ...
  binCenterDeg, pairCount, probability, probabilityDensityPerDegree, ...
  cumulativeProbability, 'VariableNames', ...
  {'bin_lower_deg','bin_upper_deg','bin_width_deg','bin_center_deg', ...
  'pair_count','probability','probability_density_per_degree', ...
  'cumulative_probability'});
```

- [ ] **Step 5: Run the synthetic tests and verify GREEN**

Run the Step 2 command.

Expected: MATLAB exits 0 and prints `NEIGHBOR_MISORIENTATION_TESTS_OK`.

- [ ] **Step 6: Commit the pure calculation unit**

```powershell
git add -- tools/mtex/four_neighbor_pairs.m tools/mtex/summarize_misorientation_angles.m tools/mtex/test_neighbor_misorientation_distribution.m
git commit -m "feat: calculate neighbor misorientation statistics"
```

### Task 2: Six-scan MTEX distribution, pooling, and figure

**Files:**
- Create: `tools/mtex/generate_neighbor_misorientation_distribution.m`
- Modify: `tools/mtex/test_neighbor_misorientation_distribution.m`
- Create: `results/mtex_neighbor_misorientation/neighbor_misorientation_distribution.csv`
- Create: `results/mtex_neighbor_misorientation/neighbor_misorientation_summary.csv`
- Create: `results/mtex_neighbor_misorientation/neighbor_misorientation_distribution.png`

**Interfaces:**
- Consumes: `scanRoot` and `outputDir` scalar strings.
- Produces: `[summary,distribution] = generate_neighbor_misorientation_distribution(scanRoot,outputDir)` and three derived files.

- [ ] **Step 1: Run the full-data test and verify integration RED**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath('C:/Users/22069/Documents/GitHub/work_hardening_paper/tools/mtex'); test_neighbor_misorientation_distribution('C:/Users/22069/Documents/GitHub/work_hardening_paper/data/ebsd_kpl_250221_7_df/scans','C:/Users/22069/Documents/GitHub/work_hardening_paper/results/mtex_neighbor_misorientation')"
```

Expected: MATLAB exits nonzero because `generate_neighbor_misorientation_distribution` is undefined.

- [ ] **Step 2: Implement six-scan neighbor-angle extraction**

Use the established folder/file/sample/deformation arrays. For each scan:

```matlab
ebsd = EBSD.load(inputFile, "convertEuler2SpatialReferenceFrame");
ebsdGrid = ebsd.gridify;
alpha = ebsdGrid("Ti-Hex");
alphaPhaseId = unique(alpha.phaseId);
assert(isscalar(alphaPhaseId));
phaseMask = reshape(ebsdGrid.phaseId, size(ebsdGrid)) == alphaPhaseId;
[index1,index2] = four_neighbor_pairs(phaseMask);
orientation1 = ebsdGrid(index1).orientations;
orientation2 = ebsdGrid(index2).orientations;
thetaDeg = double(angle(orientation1,orientation2) / degree);
[sampleStats,sampleDistribution] = ...
  summarize_misorientation_angles(thetaDeg,0:0.1:94);
```

Prepend sample metadata to each summary row and distribution block. Preserve each angle vector in a six-element cell array for pooled statistics, then call the same summary helper on `vertcat(angleCells{:})`.

- [ ] **Step 3: Validate and write CSV files**

Assert seven summary rows, 6580 distribution rows, pair-count bounds, per-sample probability conservation, and pooled bin counts equal the sum of the six sample bin counts. Write both tables with `writetable`.

- [ ] **Step 4: Render the two-panel distribution figure**

Create an invisible 1200×900 figure with `tiledlayout(2,1)`. Plot six sample probability-density curves as colored solid lines and the pooled curve as a thicker black dashed line in both panels. Set the upper x-limit to `[0,94]` and the lower x-limit to `[0,10]`; use shared labels in degrees and probability density per degree. Export at 300 dpi with a white background and remove axes toolbars before export.

- [ ] **Step 5: Run the full integration test and verify GREEN**

Run the Step 1 command.

Expected: MATLAB exits 0, prints seven summary lines and `NEIGHBOR_MISORIENTATION_TESTS_OK`, and creates both CSV files plus the PNG.

- [ ] **Step 6: Run static, structural, and visual checks**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "files={'C:/Users/22069/Documents/GitHub/work_hardening_paper/tools/mtex/four_neighbor_pairs.m','C:/Users/22069/Documents/GitHub/work_hardening_paper/tools/mtex/summarize_misorientation_angles.m','C:/Users/22069/Documents/GitHub/work_hardening_paper/tools/mtex/generate_neighbor_misorientation_distribution.m','C:/Users/22069/Documents/GitHub/work_hardening_paper/tools/mtex/test_neighbor_misorientation_distribution.m'}; for i=1:numel(files), issues=checkcode(files{i},'-id'); assert(isempty(issues),files{i}); end; disp('CHECKCODE_OK');"
git diff --check
```

Visually inspect the PNG at original resolution and verify that all seven curves are legible, the low-angle panel resolves the 0° peak, and no label, legend, or axis is clipped.

- [ ] **Step 7: Commit the verified distribution analysis**

```powershell
git add -- tools/mtex/generate_neighbor_misorientation_distribution.m tools/mtex/test_neighbor_misorientation_distribution.m results/mtex_neighbor_misorientation
git commit -m "feat: compare Ti-Hex neighbor misorientation distributions"
```

