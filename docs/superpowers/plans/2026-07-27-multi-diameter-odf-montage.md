# Multi-diameter ODF Montage Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate a publication-ready 3-by-2 montage comparing pixel-weighted alpha-Ti ODF sections for all six raw EBSD scans.

**Architecture:** A focused MATLAB entry point reuses the registered EBSD catalog and loader, calculates one normalized pixel-weighted ODF per raw scan, renders the three registered `phi2` sections, and assembles the six panels with a shared global MRD scale. A dedicated integration test validates catalog selection, numerical summaries, shared scaling, and non-empty PNG/PDF outputs.

**Tech Stack:** MATLAB R2025a, MTEX 6.1.1, project EBSD helpers, MATLAB `exportgraphics`.

## Global Constraints

- Use only the six raw CTF files registered by `comprehensive_ebsd_catalog`.
- Use indexed `Ti-Hex` pixels and pixel weighting.
- Use a De la Vallee Poussin kernel with 5 deg halfwidth.
- Render `phi2 = 0 deg`, `30 deg`, and `60 deg` at 5 deg resolution.
- Use X = AD, Y = TD/RD, and Z = ND.
- Use one global MRD range from 0 to the maximum across all six ODFs.
- Use a perceptually ordered non-rainbow color map.
- Do not modify raw or denoised CTF files.
- Export PNG, PDF, and a six-row numerical summary CSV under `results/mtex_odf_diameter_montage`.

---

### Task 1: Define the montage contract with a failing integration test

**Files:**
- Create: `tools/mtex/test_generate_odf_diameter_montage.m`
- Test: `tools/mtex/test_generate_odf_diameter_montage.m`

**Interfaces:**
- Consumes: `comprehensive_ebsd_catalog(scanRoot)`.
- Produces: the required public interface
  `summary = generate_odf_diameter_montage(scanRoot, outputRoot)`.

- [ ] **Step 1: Write the failing test**

```matlab
function test_generate_odf_diameter_montage(scanRoot, outputRoot)
arguments
  scanRoot (1,1) string
  outputRoot (1,1) string
end

assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded for this test.");
if isfolder(outputRoot)
  rmdir(outputRoot, "s");
end

summary = generate_odf_diameter_montage(scanRoot, outputRoot);
expectedColumns = ["sample", "diameter_mm", "cold_reduction_percent", ...
  "input_path", "valid_ti_hex_orientation_count", ...
  "kernel_halfwidth_deg", "grid_resolution_deg", "maximum_mrd", ...
  "global_color_limit_max_mrd"];
assert(isequal(string(summary.Properties.VariableNames), expectedColumns));
assert(height(summary) == 6);
assert(isequal(summary.sample, ...
  ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"]));
assert(isequal(summary.diameter_mm, [7;6.48;6.02;5.6;5.25;5]));
assert(isequal(summary.cold_reduction_percent, ...
  [0;14.31;26.04;36;43.75;48.98]));
assert(all(summary.valid_ti_hex_orientation_count > 0));
assert(all(summary.kernel_halfwidth_deg == 5));
assert(all(summary.grid_resolution_deg == 5));
assert(all(isfinite(summary.maximum_mrd) & summary.maximum_mrd > 0));
assert(numel(unique(summary.global_color_limit_max_mrd)) == 1);
assert(summary.global_color_limit_max_mrd(1) == max(summary.maximum_mrd));

csvPath = fullfile(outputRoot, "odf_diameter_summary.csv");
pngPath = fullfile(outputRoot, "odf_diameter_montage.png");
pdfPath = fullfile(outputRoot, "odf_diameter_montage.pdf");
assert(isfile(csvPath) && dir(csvPath).bytes > 0);
assert(isfile(pngPath) && dir(pngPath).bytes > 0);
assert(isfile(pdfPath) && dir(pdfPath).bytes > 0);
roundTrip = readtable(csvPath, "TextType", "string");
assert(isequal(string(roundTrip.Properties.VariableNames), expectedColumns));
assert(height(roundTrip) == 6);
fprintf("test_generate_odf_diameter_montage passed\n");
end
```

- [ ] **Step 2: Run the test and verify the expected RED state**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "startup_mtex('noMenu'); addpath('tools/mtex'); test_generate_odf_diameter_montage(string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')), string(fullfile(pwd,'.codex_tmp','odf-montage-red')));"
```

Expected: MATLAB exits non-zero because
`generate_odf_diameter_montage` is undefined.

- [ ] **Step 3: Commit the failing test**

```powershell
git add -- tools/mtex/test_generate_odf_diameter_montage.m
git commit -m "test: define multi-diameter ODF montage"
```

### Task 2: Implement ODF calculation, rendering, and outputs

**Files:**
- Create: `tools/mtex/generate_odf_diameter_montage.m`
- Test: `tools/mtex/test_generate_odf_diameter_montage.m`

**Interfaces:**
- Consumes:
  `catalog = comprehensive_ebsd_catalog(scanRoot)` and
  `[ebsdFull, meta] = load_comprehensive_ebsd_scan(catalogRow)`.
- Produces:
  `summary = generate_odf_diameter_montage(scanRoot, outputRoot)` and
  `odf_diameter_montage.png`, `odf_diameter_montage.pdf`,
  `odf_diameter_summary.csv`.

- [ ] **Step 1: Implement registered raw-scan selection and ODF calculation**

Create `generate_odf_diameter_montage.m` with an arguments block, MTEX/input
assertions, output-directory creation, and:

```matlab
catalog = comprehensive_ebsd_catalog(scanRoot);
catalog = catalog(catalog.variant == "raw", :);
assert(height(catalog) == 6);
assert(isequal(catalog.sample, ...
  ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"]));

kernelHalfwidthDeg = 5;
gridResolutionDeg = 5;
kernel = SO3DeLaValleePoussinKernel( ...
  "halfwidth", kernelHalfwidthDeg * degree);
odfModels = cell(6,1);
orientationCounts = zeros(6,1);
maximumMrd = zeros(6,1);

for scanIndex = 1:6
  [ebsdFull, ~] = load_comprehensive_ebsd_scan(catalog(scanIndex,:));
  tiEbsd = ebsdFull("Ti-Hex");
  assert(~isempty(tiEbsd));
  orientationCounts(scanIndex) = length(tiEbsd);
  rbfOdf = calcDensity(tiEbsd.orientations, "kernel", kernel, ...
    "weights", ones(length(tiEbsd),1), "silent");
  normalizedRbfOdf = rbfOdf / double(mean(rbfOdf));
  odfModels{scanIndex} = SO3FunHarmonic(normalizedRbfOdf, ...
    "bandwidth", kernel.bandwidth);
  maximumMrd(scanIndex) = double(max(odfModels{scanIndex}, ...
    "resolution", gridResolutionDeg * degree));
end
globalMaximumMrd = max(maximumMrd);
```

Assert that means are one within `1e-6`, counts are positive, and maxima are
finite and positive.

- [ ] **Step 2: Implement panel rendering and 3-by-2 montage assembly**

Add private local helpers that:

1. render each ODF into a temporary RGB panel with:

```matlab
plotSection(odf, "phi2", [0 30 60] * degree, ...
  "contourf", "silent", "layout", [1 3], ...
  "resolution", gridResolutionDeg * degree, ...
  "colorRange", [0 globalMaximumMrd]);
mtexColorMap parula
```

2. crop white margins without changing the plotted content;
3. place the six RGB panels into `tiledlayout(3,2)` in registered order;
4. title each tile with diameter and cold reduction;
5. attach one east-side colorbar labelled `ODF intensity (MRD)` with
   limits `[0 globalMaximumMrd]`;
6. export:

```matlab
exportgraphics(figureHandle, pngPath, "Resolution", 600, ...
  "BackgroundColor", "white");
exportgraphics(figureHandle, pdfPath, "ContentType", "image", ...
  "BackgroundColor", "white");
```

Temporary files must be owned by a `tempname` directory and removed through
an `onCleanup` handler.

- [ ] **Step 3: Write the six-row numerical summary**

Construct the table with exact schema:

```matlab
summary = table(catalog.sample, catalog.diameter_mm, ...
  catalog.cold_reduction_percent, catalog.input_path, ...
  orientationCounts, repmat(kernelHalfwidthDeg,6,1), ...
  repmat(gridResolutionDeg,6,1), maximumMrd, ...
  repmat(globalMaximumMrd,6,1), ...
  "VariableNames", ["sample", "diameter_mm", ...
  "cold_reduction_percent", "input_path", ...
  "valid_ti_hex_orientation_count", "kernel_halfwidth_deg", ...
  "grid_resolution_deg", "maximum_mrd", ...
  "global_color_limit_max_mrd"]);
writetable(summary, fullfile(outputRoot, "odf_diameter_summary.csv"));
```

Assert all three output files exist and are non-empty before returning.

- [ ] **Step 4: Run the integration test and verify GREEN**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "startup_mtex('noMenu'); addpath('tools/mtex'); test_generate_odf_diameter_montage(string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')), string(fullfile(pwd,'.codex_tmp','odf-montage-test')));"
```

Expected: exit code 0 and
`test_generate_odf_diameter_montage passed`.

- [ ] **Step 5: Commit the implementation**

```powershell
git add -- tools/mtex/generate_odf_diameter_montage.m
git commit -m "feat: generate multi-diameter ODF montage"
```

### Task 3: Generate the final figure and verify manuscript readiness

**Files:**
- Create:
  `results/mtex_odf_diameter_montage/odf_diameter_montage.png`
- Create:
  `results/mtex_odf_diameter_montage/odf_diameter_montage.pdf`
- Create:
  `results/mtex_odf_diameter_montage/odf_diameter_summary.csv`

**Interfaces:**
- Consumes:
  `generate_odf_diameter_montage(scanRoot, outputRoot)`.
- Produces: the three final derived artifacts under the registered results
  directory.

- [ ] **Step 1: Run the final generator from a fresh MATLAB process**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "startup_mtex('noMenu'); addpath('tools/mtex'); generate_odf_diameter_montage(string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')), string(fullfile(pwd,'results','mtex_odf_diameter_montage'));"
```

Expected: exit code 0 and three non-empty final artifacts.

- [ ] **Step 2: Verify numerical and file-level acceptance criteria**

```powershell
Get-Item 'results\mtex_odf_diameter_montage\odf_diameter_montage.png',
  'results\mtex_odf_diameter_montage\odf_diameter_montage.pdf',
  'results\mtex_odf_diameter_montage\odf_diameter_summary.csv' |
  Select-Object FullName,Length
Import-Csv 'results\mtex_odf_diameter_montage\odf_diameter_summary.csv' |
  Format-Table sample,diameter_mm,cold_reduction_percent,
    valid_ti_hex_orientation_count,maximum_mrd,
    global_color_limit_max_mrd
```

Expected: six rows, positive file sizes/counts/maxima, and one identical
global color-limit value.

- [ ] **Step 3: Visually inspect the PNG**

Open the PNG and verify six ordered diameter panels, three complete `phi2`
sections per panel, readable labels, one shared MRD colorbar, no blank or
clipped plots, and consistent color limits.

- [ ] **Step 4: Commit final derived artifacts**

```powershell
git add -- results/mtex_odf_diameter_montage
git commit -m "results: add multi-diameter ODF montage"
```
