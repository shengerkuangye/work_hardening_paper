# α-Ti C-Axis AD Alignment Metrics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute, organize, and visualize area-weighted α-Ti c-axis alignment relative to the bar axial direction for all six EBSD scans.

**Architecture:** A pure MATLAB statistics function consumes acute c-axis-to-AD angles and produces validated scalar metrics plus a 2° binned distribution. A separate MTEX batch function imports the six original CTF files, converts `[0001]` directions to specimen coordinates, calls the statistics function, writes two CSV files, and renders two publication-ready comparison figures.

**Tech Stack:** MATLAB R2025a, MTEX 6.1.1, MATLAB tables/graphics, CTF EBSD input.

## Global Constraints

- Do not modify or overwrite any raw CTF data.
- Use the six original CTF files selected by the existing pole-figure workflow; exclude `_denoised.ctf` files.
- Use only indexed `Ti-Hex` pixels and equal per-pixel weights.
- AD is specimen x, radial direction is specimen y, and scan-plane normal is specimen z.
- Compute the antipodal acute angle in the closed interval 0°–90°.
- Use fixed 2° bins from 0° through 90° for every sample.
- Write all derived files under `results/mtex_c_axis_ad_alignment/`.
- Preserve the sample order `7d`, `6.48d`, `6.02d`, `5.6d`, `5.25d`, `5d`.

---

## File Structure

- Create `tools/mtex/compute_c_axis_ad_statistics.m`: pure numerical statistics and histogram function; no file I/O.
- Create `tools/mtex/generate_c_axis_ad_alignment_metrics.m`: MTEX import, c-axis conversion, CSV assembly, plotting, and export.
- Create `tools/mtex/test_c_axis_ad_alignment_metrics.m`: synthetic unit checks and optional full-data integration checks.
- Create `results/mtex_c_axis_ad_alignment/c_axis_ad_alignment_summary.csv`: six-row scalar summary.
- Create `results/mtex_c_axis_ad_alignment/c_axis_ad_angle_distribution.csv`: long-form 270-row distribution table.
- Create `results/mtex_c_axis_ad_alignment/c_axis_ad_alignment_metrics.png`: three-panel scalar metric comparison.
- Create `results/mtex_c_axis_ad_alignment/c_axis_ad_angle_distribution.png`: six measured distributions plus random reference.

### Task 1: Pure numerical alignment statistics

**Files:**
- Create: `tools/mtex/compute_c_axis_ad_statistics.m`
- Create: `tools/mtex/test_c_axis_ad_alignment_metrics.m`

**Interfaces:**
- Consumes: `thetaDeg` as a finite numeric vector in `[0,90]`, `totalPixelCount` as a positive integer not less than `numel(thetaDeg)`, and `binEdgesDeg` as increasing numeric edges spanning `[0,90]`.
- Produces: `[stats, distribution] = compute_c_axis_ad_statistics(thetaDeg, totalPixelCount, binEdgesDeg)`, where `stats` is a scalar struct and `distribution` is a table.

- [ ] **Step 1: Write the failing synthetic test**

Create `tools/mtex/test_c_axis_ad_alignment_metrics.m` with:

```matlab
function test_c_axis_ad_alignment_metrics
edges = 0:2:90;

[parallelStats, parallelDist] = compute_c_axis_ad_statistics( ...
  zeros(100,1), 100, edges);
assert(abs(parallelStats.alignment_factor - 1) < 1e-12);
assert(parallelStats.fraction_within_15deg == 1);
assert(sum(parallelDist.pixel_count) == 100);
assert(abs(sum(parallelDist.probability) - 1) < 1e-12);

[normalStats, normalDist] = compute_c_axis_ad_statistics( ...
  90 * ones(100,1), 100, edges);
assert(abs(normalStats.alignment_factor + 0.5) < 1e-12);
assert(normalStats.fraction_within_45deg == 0);
assert(sum(normalDist.pixel_count) == 100);

[mixedStats, ~] = compute_c_axis_ad_statistics([0; 90], 4, edges);
assert(abs(mixedStats.mean_cos2 - 0.5) < 1e-12);
assert(abs(mixedStats.alignment_factor - 0.25) < 1e-12);
assert(abs(mixedStats.alpha_ti_pixel_fraction - 0.5) < 1e-12);

fprintf("C_AXIS_AD_ALIGNMENT_TESTS_OK\n");
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath('C:/Users/22069/Documents/GitHub/work_hardening_paper/tools/mtex'); test_c_axis_ad_alignment_metrics"
```

Expected: MATLAB exits nonzero because `compute_c_axis_ad_statistics` is undefined.

- [ ] **Step 3: Implement the pure statistics function**

Create `tools/mtex/compute_c_axis_ad_statistics.m` with input validation, population statistics, fixed-threshold proportions, and `histcounts` output. The calculation body must use:

```matlab
thetaDeg = double(thetaDeg(:));
n = numel(thetaDeg);
assert(n > 0, "thetaDeg must not be empty.");
assert(all(isfinite(thetaDeg)), "thetaDeg contains nonfinite values.");
assert(all(thetaDeg >= 0 & thetaDeg <= 90), ...
  "thetaDeg values must be between 0 and 90 degrees.");
assert(totalPixelCount >= n && totalPixelCount == fix(totalPixelCount), ...
  "totalPixelCount must be an integer not less than the valid pixel count.");
assert(binEdgesDeg(1) == 0 && binEdgesDeg(end) == 90 && ...
  all(diff(binEdgesDeg) > 0), ...
  "binEdgesDeg must increase from 0 to 90 degrees.");

cos2 = cosd(thetaDeg).^2;
stats.total_pixel_count = totalPixelCount;
stats.alpha_ti_pixel_count = n;
stats.alpha_ti_pixel_fraction = n / totalPixelCount;
stats.mean_angle_deg = mean(thetaDeg);
stats.std_angle_deg = std(thetaDeg, 1);
quartiles = prctile(thetaDeg, [25, 50, 75]);
stats.q25_angle_deg = quartiles(1);
stats.median_angle_deg = quartiles(2);
stats.q75_angle_deg = quartiles(3);
stats.mean_cos2 = mean(cos2);
stats.alignment_factor = (3 * stats.mean_cos2 - 1) / 2;
stats.fraction_within_15deg = mean(thetaDeg <= 15);
stats.fraction_within_30deg = mean(thetaDeg <= 30);
stats.fraction_within_45deg = mean(thetaDeg <= 45);

pixelCount = histcounts(thetaDeg, binEdgesDeg)';
binLowerDeg = binEdgesDeg(1:end-1)';
binUpperDeg = binEdgesDeg(2:end)';
binWidthDeg = binUpperDeg - binLowerDeg;
binCenterDeg = (binLowerDeg + binUpperDeg) / 2;
probability = pixelCount / n;
probabilityDensityPerDegree = probability ./ binWidthDeg;
distribution = table(binLowerDeg, binUpperDeg, binWidthDeg, ...
  binCenterDeg, pixelCount, probability, probabilityDensityPerDegree, ...
  'VariableNames', {'bin_lower_deg','bin_upper_deg','bin_width_deg', ...
  'bin_center_deg','pixel_count','probability', ...
  'probability_density_per_degree'});
```

- [ ] **Step 4: Run the synthetic test to verify it passes**

Run the command from Step 2.

Expected: MATLAB exits 0 and prints `C_AXIS_AD_ALIGNMENT_TESTS_OK`.

- [ ] **Step 5: Commit the numerical unit**

```powershell
git add -- tools/mtex/compute_c_axis_ad_statistics.m tools/mtex/test_c_axis_ad_alignment_metrics.m
git commit -m "feat: calculate c-axis axial alignment statistics"
```

### Task 2: Six-scan MTEX batch analysis and CSV outputs

**Files:**
- Create: `tools/mtex/generate_c_axis_ad_alignment_metrics.m`
- Modify: `tools/mtex/test_c_axis_ad_alignment_metrics.m`
- Create: `results/mtex_c_axis_ad_alignment/c_axis_ad_alignment_summary.csv`
- Create: `results/mtex_c_axis_ad_alignment/c_axis_ad_angle_distribution.csv`

**Interfaces:**
- Consumes: `scanRoot` and `outputDir` scalar strings.
- Produces: `[summary, distribution] = generate_c_axis_ad_alignment_metrics(scanRoot, outputDir)` and two CSV files.

- [ ] **Step 1: Extend the test with CSV integration assertions**

Replace the Task 1 function declaration with the following declaration and arguments block, leaving the existing synthetic assertions immediately after it:

```matlab
function test_c_axis_ad_alignment_metrics(scanRoot, outputDir)
arguments
  scanRoot (1,1) string = ""
  outputDir (1,1) string = ""
end
```

Then insert this block immediately before `fprintf("C_AXIS_AD_ALIGNMENT_TESTS_OK\n");`:

```matlab
if scanRoot ~= ""
  assert(outputDir ~= "", "outputDir is required for integration testing.");
  [summary, distribution] = generate_c_axis_ad_alignment_metrics( ...
    scanRoot, outputDir);
  assert(height(summary) == 6);
  assert(height(distribution) == 270);
  assert(all(summary.total_pixel_count == 360000));
  assert(all(summary.alpha_ti_pixel_count > 0));
  assert(all(summary.alignment_factor >= -0.5 & ...
    summary.alignment_factor <= 1));
  assert(all(summary.fraction_within_15deg >= 0 & ...
    summary.fraction_within_15deg <= 1));
  for sample = unique(distribution.sample, "stable")'
    rows = distribution.sample == sample;
    assert(abs(sum(distribution.probability(rows)) - 1) < 1e-10);
    assert(abs(sum(distribution.probability_density_per_degree(rows) .* ...
      distribution.bin_width_deg(rows)) - 1) < 1e-10);
  end
  assert(isfile(fullfile(outputDir, "c_axis_ad_alignment_summary.csv")));
  assert(isfile(fullfile(outputDir, "c_axis_ad_angle_distribution.csv")));
end
```

- [ ] **Step 2: Run the integration test to verify the batch function is missing**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath('C:/Users/22069/Documents/GitHub/work_hardening_paper/tools/mtex'); test_c_axis_ad_alignment_metrics('C:/Users/22069/Documents/GitHub/work_hardening_paper/data/ebsd_kpl_250221_7_df/scans','C:/Users/22069/Documents/GitHub/work_hardening_paper/results/mtex_c_axis_ad_alignment')"
```

Expected: MATLAB exits nonzero because `generate_c_axis_ad_alignment_metrics` is undefined.

- [ ] **Step 3: Implement fixed sample mapping and c-axis angle conversion**

The batch function must use the same folder/file/sample/deformation arrays as `generate_c_axis_pole_figures.m`, load each CTF with `convertEuler2SpatialReferenceFrame`, select `Ti-Hex`, and compute:

```matlab
cAxis = Miller(0, 0, 0, 1, alpha.CS);
cDirections = alpha.orientations * cAxis;
cosToAD = min(1, abs(cDirections.x));
thetaDeg = acosd(cosToAD);
[sampleStats, sampleDistribution] = compute_c_axis_ad_statistics( ...
  thetaDeg, length(ebsd), 0:2:90);
```

For each sample, prepend `sample`, `folder`, `input_file`, and `cold_reduction_percent` to the scalar summary row and prepend `sample`, `cold_reduction_percent` to all 45 distribution rows. Add repeated theoretical-reference columns to the summary:

```matlab
randomMeanAngleDeg = 180 / pi;
randomMedianAngleDeg = 60;
randomAlignmentFactor = 0;
randomFraction15 = 1 - cosd(15);
randomFraction30 = 1 - cosd(30);
randomFraction45 = 1 - cosd(45);
```

- [ ] **Step 4: Write the two CSV files**

Use:

```matlab
writetable(summary, fullfile(outputDir, ...
  "c_axis_ad_alignment_summary.csv"));
writetable(distribution, fullfile(outputDir, ...
  "c_axis_ad_angle_distribution.csv"));
```

Assert six summary rows, 270 distribution rows, finite metrics, valid ranges, per-sample bin counts equal Ti-Hex pixel counts, and per-sample probabilities sum to one.

- [ ] **Step 5: Run the CSV integration test to verify it passes**

Run the integration command from Step 2.

Expected: MATLAB exits 0, creates both CSV files, and prints `C_AXIS_AD_ALIGNMENT_TESTS_OK`.

- [ ] **Step 6: Commit the verified CSV analysis**

```powershell
git add -- tools/mtex/generate_c_axis_ad_alignment_metrics.m results/mtex_c_axis_ad_alignment/c_axis_ad_alignment_summary.csv results/mtex_c_axis_ad_alignment/c_axis_ad_angle_distribution.csv
git commit -m "feat: summarize c-axis axial alignment metrics"
```

### Task 3: Comparative figures and full verification

**Files:**
- Modify: `tools/mtex/generate_c_axis_ad_alignment_metrics.m`
- Create: `results/mtex_c_axis_ad_alignment/c_axis_ad_alignment_metrics.png`
- Create: `results/mtex_c_axis_ad_alignment/c_axis_ad_angle_distribution.png`

**Interfaces:**
- Consumes: `summary` and `distribution` tables assembled by Task 2.
- Produces: two 300-dpi PNG figures with fixed dimensions and no interactive axes toolbar.

- [ ] **Step 1: Add failing figure-output assertions**

Inside the integration block of `test_c_axis_ad_alignment_metrics.m`, add:

```matlab
assert(isfile(fullfile(outputDir, "c_axis_ad_alignment_metrics.png")));
assert(isfile(fullfile(outputDir, "c_axis_ad_angle_distribution.png")));
```

Run the full integration command. Expected: MATLAB exits nonzero because `c_axis_ad_alignment_metrics.png` does not exist.

- [ ] **Step 2: Implement the three-panel scalar comparison**

Create an invisible 1200×1000 figure and use `tiledlayout(3,1)`. Plot:

```matlab
x = summary.cold_reduction_percent;

% Panel a
plot(x, summary.alignment_factor, "-o", "LineWidth", 1.5);
yline(0, "--", "Random");
ylabel("F_{AD}");

% Panel b
plot(x, summary.mean_angle_deg, "-o", "LineWidth", 1.5);
hold on;
plot(x, summary.median_angle_deg, "-s", "LineWidth", 1.5);
yline(180/pi, "--", "Random mean");
yline(60, ":", "Random median");
ylabel("Angle to AD (degree)");

% Panel c
plot(x, 100 * summary.fraction_within_15deg, "-o", ...
  x, 100 * summary.fraction_within_30deg, "-s", ...
  x, 100 * summary.fraction_within_45deg, "-^", ...
  "LineWidth", 1.5);
ylabel("Fraction (%)");
xlabel("Cold reduction (%)");
```

Use sample labels as data labels, add legends, grids, panel labels `(a)`–`(c)`, and export at 300 dpi with a white background.

- [ ] **Step 3: Implement the angle-distribution comparison**

Plot the six `probability_density_per_degree` series against `bin_center_deg`. Add the theoretical random curve:

```matlab
thetaReference = linspace(0, 90, 361)';
randomDensityPerDegree = (pi/180) * sind(thetaReference);
plot(thetaReference, randomDensityPerDegree, "k--", ...
  "LineWidth", 1.8, "DisplayName", "Random");
```

Fix x limits to `[0,90]`, label axes `c-axis angle to AD (degree)` and `Probability density (degree^{-1})`, add a legend and grid, then export at 300 dpi with a white background.

- [ ] **Step 4: Run the full integration test**

Run the Task 2 Step 1 command.

Expected: MATLAB exits 0 and prints `C_AXIS_AD_ALIGNMENT_TESTS_OK` after importing all six scans, writing both CSVs, and exporting both figures.

- [ ] **Step 5: Run static and repository checks**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "files={'C:/Users/22069/Documents/GitHub/work_hardening_paper/tools/mtex/compute_c_axis_ad_statistics.m','C:/Users/22069/Documents/GitHub/work_hardening_paper/tools/mtex/generate_c_axis_ad_alignment_metrics.m','C:/Users/22069/Documents/GitHub/work_hardening_paper/tools/mtex/test_c_axis_ad_alignment_metrics.m'}; for i=1:numel(files), issues=checkcode(files{i},'-id'); assert(isempty(issues),files{i}); end; disp('CHECKCODE_OK');"
git diff --check
```

Expected: MATLAB prints `CHECKCODE_OK`; `git diff --check` exits 0.

- [ ] **Step 6: Visually inspect both PNG files**

Verify all sample labels are legible, the scalar figure contains three panels, the distribution figure contains six sample curves plus the random reference, axes are not clipped, and no axes toolbar appears.

- [ ] **Step 7: Commit the figures and final generator**

```powershell
git add -- tools/mtex/generate_c_axis_ad_alignment_metrics.m tools/mtex/test_c_axis_ad_alignment_metrics.m results/mtex_c_axis_ad_alignment
git commit -m "feat: plot c-axis axial alignment evolution"
```
