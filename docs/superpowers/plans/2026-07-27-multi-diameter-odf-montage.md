# Multi-diameter ODF Diagnostic Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate a 6-row by 7-column alpha-Ti ODF diagnostic for six raw EBSD scans, with six-row sample and 42-row sample--section MRD summaries for later manuscript-section selection.

**Architecture:** A MATLAB entry point reuses the registered catalog and loader, estimates one pixel-weighted ODF per scan with `6/mmm` and `SS = 1`, evaluates `phi2 = 0:10:60 deg`, writes both summaries, derives one global MRD limit from all 42 section maxima, then renders six identical seven-section row rasters and composes them at exact pixel size. It does not select manuscript representative sections.

**Tech Stack:** MATLAB R2025a, MTEX 6.1.1, project EBSD helpers, MATLAB `exportgraphics`.

## Global Constraints

- Use only six raw CTF files from `comprehensive_ebsd_catalog` and indexed `Ti-Hex` pixels with pixel weighting.
- Use `6/mmm` crystal symmetry, the `622` rotational fundamental zone, and bar-specimen `SS = 1`; never apply `222`.
- Use one 5 deg De la Vallee Poussin kernel, 5 deg evaluation resolution, and `phi2 = 0:10:60 deg` for every sample.
- Use one MRD range from zero to the maximum across all 42 section maxima and a perceptually ordered non-rainbow map.
- Do not modify raw or denoised CTF files.
- Export `odf_diameter_full_sections.png`, `odf_diameter_full_sections.pdf`, `odf_diameter_summary.csv` (six rows), and `odf_diameter_section_summary.csv` (42 rows) under `results/mtex_odf_diameter_montage`.
- Preserve approximately 600 dpi in the PNG metadata and at least 590 dpi of page-equivalent raster data in the image-content PDF.
- Select manuscript angles only after checking the full diagnostic and the 42-row CSV. A later caption must use `selected phi2 sections` and list their angles.

---

### Task 1: Define the output contract with a failing integration test

**Files:**

- Create: `tools/mtex/test_generate_odf_diameter_montage.m`
- Test: `tools/mtex/test_generate_odf_diameter_montage.m`

**Interface:** `[sampleSummary, sectionSummary] = generate_odf_diameter_montage(scanRoot, outputRoot)`.

- [ ] **Step 1: Write the failing test**

```matlab
function test_generate_odf_diameter_montage(scanRoot)
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded for this test.");
outputRoot = string(tempname);
mkdir(outputRoot);
cleanupOutput = onCleanup(@() rmdir(outputRoot,"s"));
[sampleSummary, sectionSummary] = generate_odf_diameter_montage(scanRoot, outputRoot);
sampleColumns = ["sample","diameter_mm","cold_reduction_percent", ...
  "input_path","valid_ti_hex_orientation_count","crystal_symmetry", ...
  "rotational_fundamental_zone","specimen_symmetry", ...
  "kernel_halfwidth_deg","grid_resolution_deg","maximum_section_mrd", ...
  "global_color_limit_max_mrd"];
sectionColumns = ["sample","diameter_mm","cold_reduction_percent", ...
  "input_path","valid_ti_hex_orientation_count","crystal_symmetry", ...
  "rotational_fundamental_zone","specimen_symmetry", ...
  "kernel_halfwidth_deg","grid_resolution_deg","phi2_deg", ...
  "section_maximum_mrd","global_color_limit_max_mrd"];
assert(isequal(string(sampleSummary.Properties.VariableNames),sampleColumns));
assert(isequal(string(sectionSummary.Properties.VariableNames),sectionColumns));
assert(height(sampleSummary) == 6);
assert(height(sectionSummary) == 42);
assert(isequal(sampleSummary.sample,["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"]));
assert(isequal(sampleSummary.diameter_mm,[7;6.48;6.02;5.6;5.25;5]));
for sampleIndex = 1:6
  sectionRows = (sampleIndex - 1) * 7 + (1:7);
  assert(isequal(sectionSummary.phi2_deg(sectionRows),(0:10:60)'));
end
assert(all(sampleSummary.crystal_symmetry == "6/mmm"));
assert(all(sampleSummary.rotational_fundamental_zone == "622"));
assert(all(sampleSummary.specimen_symmetry == "1"));
assert(all(sampleSummary.kernel_halfwidth_deg == 5));
assert(all(sampleSummary.grid_resolution_deg == 5));
assert(all(isfinite(sectionSummary.section_maximum_mrd) & sectionSummary.section_maximum_mrd > 0));
assert(numel(unique(sectionSummary.global_color_limit_max_mrd)) == 1);
assert(sectionSummary.global_color_limit_max_mrd(1) == max(sectionSummary.section_maximum_mrd));
assert(isfile(fullfile(outputRoot,"odf_diameter_summary.csv")));
assert(isfile(fullfile(outputRoot,"odf_diameter_section_summary.csv")));
assert(isfile(fullfile(outputRoot,"odf_diameter_full_sections.png")));
assert(isfile(fullfile(outputRoot,"odf_diameter_full_sections.pdf")));
roundTripSample = readtable(fullfile(outputRoot,"odf_diameter_summary.csv"),"TextType","string");
roundTripSection = readtable(fullfile(outputRoot,"odf_diameter_section_summary.csv"),"TextType","string");
assert(isequal(string(roundTripSample.Properties.VariableNames),sampleColumns));
assert(isequal(string(roundTripSection.Properties.VariableNames),sectionColumns));
assert(height(roundTripSample) == 6 && height(roundTripSection) == 42);
end
```

- [ ] **Step 2: Run RED**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "startup_mtex('noMenu'); addpath('tools/mtex'); test_generate_odf_diameter_montage(string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')));"
```

Expected: MATLAB exits non-zero because `generate_odf_diameter_montage` is undefined.

- [ ] **Step 3: Commit the failing test**

```powershell
git add -- tools/mtex/test_generate_odf_diameter_montage.m
git commit -m "test: define seven-section ODF diagnostic"
```

### Task 2: Implement ODF evaluation and both CSV outputs

**Files:**

- Create: `tools/mtex/generate_odf_diameter_montage.m`
- Create: `tools/mtex/normalize_positive_mean_density.m`
- Test: `tools/mtex/test_generate_odf_diameter_montage.m`

**Interface:** Returns `sampleSummary` (six rows) and `sectionSummary` (42 rows), and writes their CSV counterparts.

- [ ] **Step 1: Implement catalog selection and validate the imported orientation symmetries**

Set `diameterOrder = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"]`, `phi2Deg = (0:10:60)'`, `kernelHalfwidthDeg = 5`, and `gridResolutionDeg = 5`. Select only raw rows from `comprehensive_ebsd_catalog(scanRoot)` and assert six rows in `diameterOrder`. For each loaded `Ti-Hex` EBSD subset, preserve its imported orientation object and use it as the sole source of symmetry:

```matlab
tiOrientations = tiEbsd.orientations;
assert(~isempty(tiOrientations));
assert(string(tiOrientations.CS.pointGroup) == "6/mmm");
assert(string(tiOrientations.CS.properGroup.pointGroup) == "622");
assert(string(tiOrientations.SS.pointGroup) == "1");
eulerRegion = fundamentalRegionEuler(tiOrientations.CS,tiOrientations.SS);
assert(isequal(eulerRegion,[360 90 60] * degree));
```

Do not instantiate a separate `crystalSymmetry("6/mmm")` or `specimenSymmetry("1")` object, and do not assign a new symmetry to imported orientations. `fundamentalRegionEuler` is the local helper that returns the `[phi1Max PhiMax phi2Max]` Euler-region limits; it must verify `360/90/60 deg` for the imported `6/mmm`/`SS = 1` orientation pair and reject any `222` specimen symmetry.

```matlab
function eulerLimits = fundamentalRegionEuler(cs,ss)
assert(string(cs.pointGroup) == "6/mmm");
assert(string(cs.properGroup.pointGroup) == "622");
assert(string(ss.pointGroup) == "1");
eulerLimits = [360 90 60] * degree;
end
```

- [ ] **Step 2: Implement one ODF and an executable seven-section maximum per sample**

For each catalog row, calculate pixel-weighted density directly from `tiOrientations` with `SO3DeLaValleePoussinKernel("halfwidth",5*degree)`. Normalize the radial-basis and converted harmonic ODFs with `normalize_positive_mean_density`, which rejects non-scalar, complex, non-finite, or non-positive means and verifies the output mean is unity. Assert the converted ODF's actual `6/mmm`, `622`, and `SS = 1` symmetries before retaining it. For every `phi2Deg`, compute the section maximum by literal Euler-grid sampling:

```matlab
phi1Deg = 0:5:355;
PhiDeg = 0:5:90;
[phi1Grid,PhiGrid] = meshgrid(phi1Deg,PhiDeg);
sectionOrientations = orientation.byEuler(phi1Grid(:)*degree, ...
  PhiGrid(:)*degree,phi2Deg(sectionIndex)*degree, ...
  tiOrientations.CS,tiOrientations.SS);
sectionMrd = real(eval(odf,sectionOrientations));
sectionMrd = sectionMrd(:);
assert(~isempty(sectionMrd) && all(isfinite(sectionMrd)));
sectionMaximumMrd(scanIndex,sectionIndex) = max(sectionMrd);
```

Use the same `phi1 = 0:5:355`, `Phi = 0:5:90`, and `phi2 = 0:10:60 deg` grids for every sample. Assert every recorded maximum is finite and positive, then set `globalMaximumMrd = max(sectionMaximumMrd,[],"all")`.

- [ ] **Step 3: Implement exact table schemas and CSV writes**

`sampleSummary` has exactly these 12 columns, in this order: `sample`, `diameter_mm`, `cold_reduction_percent`, `input_path`, `valid_ti_hex_orientation_count`, `crystal_symmetry`, `rotational_fundamental_zone`, `specimen_symmetry`, `kernel_halfwidth_deg`, `grid_resolution_deg`, `maximum_section_mrd`, and `global_color_limit_max_mrd`.

`sectionSummary` has exactly these 13 columns, in this order: `sample`, `diameter_mm`, `cold_reduction_percent`, `input_path`, `valid_ti_hex_orientation_count`, `crystal_symmetry`, `rotational_fundamental_zone`, `specimen_symmetry`, `kernel_halfwidth_deg`, `grid_resolution_deg`, `phi2_deg`, `section_maximum_mrd`, and `global_color_limit_max_mrd`. Build the 42 rows using `repelem(catalog.sample,7)`, `repmat(phi2Deg,6,1)`, and `reshape(sectionMaximumMrd.',[],1)`, so records are ordered by diameter then phi2. Write the tables with:

```matlab
writetable(sampleSummary,fullfile(outputRoot,"odf_diameter_summary.csv"));
writetable(sectionSummary,fullfile(outputRoot,"odf_diameter_section_summary.csv"));
```

Assert six and 42 rows, the phi2 sequence for each sample, and `globalMaximumMrd == max(sectionSummary.section_maximum_mrd)`.

- [ ] **Step 4: Run the test and verify calculation assertions pass**

Run the Task 1 command using `.codex_tmp/odf-montage-summary-green`. Expected: only the still-missing PNG/PDF assertions fail; symmetry, ordering, MRD, and both CSV assertions pass.

- [ ] **Step 5: Commit calculation and summaries**

```powershell
git add -- tools/mtex/generate_odf_diameter_montage.m tools/mtex/test_generate_odf_diameter_montage.m
git commit -m "feat: calculate seven-section alpha-Ti ODF diagnostics"
```

### Task 3: Render the 6-row by 7-column diagnostic

**Files:**

- Modify: `tools/mtex/generate_odf_diameter_montage.m`
- Test: `tools/mtex/test_generate_odf_diameter_montage.m`

**Interface:** Writes the diagnostic PNG and PDF required by Task 1.

- [ ] **Step 1: Extend the test, then run RED**

Add these checks, then run the Task 1 command using `.codex_tmp/odf-montage-render-red`:

```matlab
imageInfo = imfinfo(fullfile(outputRoot,"odf_diameter_full_sections.png"));
assert(imageInfo.Width > 0 && imageInfo.Height > 0);
assert(imageInfo.Width > imageInfo.Height, ...
  "The 6-row by 7-column SS=1 diagnostic must be horizontally oriented.");
assert(imageInfo.Width >= 2400 && imageInfo.Height >= 1800, ...
  "The diagnostic raster is smaller than the planned publication output.");
```

Expected: failure because the renderer has not written the PNG/PDF.

- [ ] **Step 2: Implement the matrix renderer**

For each diameter, render the seven requested sections as one identically sized temporary row raster using `plotSection(...,"layout",[1 7])`. Require seven non-empty MTEX `sphericalPlot` axes per row and 42 rendered axes overall. Export every row at 600 dpi, crop all six rows to one common non-white bounding box, then compose the six rasters at their exact pixel size in registered diameter order:

```matlab
plotSection(odfModels{scanIndex},"phi2",phi2Deg*degree, ...
  "contourf","silent","layout",[1 7], ...
  "resolution",gridResolutionDeg*degree, ...
  "colorRange",[0 globalMaximumMrd]);
mtexColorMap parula
renderedSectionCount = renderedSectionCount + 7;
```

After the row loop, require `assert(renderedSectionCount == 42)`. Preserve the imported-symmetry assertions (`6/mmm`, proper group `622`, and `SS = 1`), add diameter/cold-reduction labels at left and `phi2` headings above the seven columns, then attach one east-side `ODF intensity (MRD)` colorbar with `[0 globalMaximumMrd]`. Export the PNG at 600 dpi. Export the PDF with `ContentType="image"` and `Resolution=600`; the test must reject PDF page-raster density below 590 dpi. Assert both outputs are non-empty.

- [ ] **Step 3: Run GREEN and commit**

Run the Task 1 command using `.codex_tmp/odf-montage-test`. Expected: exit code 0. Then run:

```powershell
git add -- tools/mtex/generate_odf_diameter_montage.m tools/mtex/test_generate_odf_diameter_montage.m
git commit -m "feat: render multi-diameter ODF diagnostic matrix"
```

### Task 4: Generate, inspect, and gate manuscript-section selection

**Files:**

- Create: `results/mtex_odf_diameter_montage/odf_diameter_full_sections.png`
- Create: `results/mtex_odf_diameter_montage/odf_diameter_full_sections.pdf`
- Create: `results/mtex_odf_diameter_montage/odf_diameter_summary.csv`
- Create: `results/mtex_odf_diameter_montage/odf_diameter_section_summary.csv`

- [ ] **Step 1: Generate final diagnostics**

Run in a fresh MATLAB process:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "startup_mtex('noMenu'); addpath('tools/mtex'); generate_odf_diameter_montage(string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')), string(fullfile(pwd,'results','mtex_odf_diameter_montage')));"
```

Expected: MATLAB exits 0 and writes four non-empty artifacts.

- [ ] **Step 2: Verify acceptance criteria**

Use `Import-Csv` to assert six sample rows, 42 section rows, unique values `6/mmm`, `622`, and `1` for the three symmetry columns, and exactly one global color-limit value equal to the maximum `section_maximum_mrd`.

- [ ] **Step 3: Inspect before selecting manuscript angles**

Open `odf_diameter_full_sections.png` and verify six ordered diameter rows, seven `phi2 = 0:10:60 deg` columns, readable labels, one shared colorbar, and no blank or clipped cells. Compare visible peaks with the 42-row summary before deciding any manuscript subset. The caption of any later manuscript figure must call the result `selected phi2 sections` and list the angles.

- [ ] **Step 4: Commit final artifacts**

```powershell
git add -- results/mtex_odf_diameter_montage
git commit -m "results: add multi-diameter ODF diagnostic"
```
