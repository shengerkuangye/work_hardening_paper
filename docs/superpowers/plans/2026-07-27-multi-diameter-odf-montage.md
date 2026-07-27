# Multi-diameter ODF Diagnostic Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate a 6-row by 7-column alpha-Ti ODF diagnostic for six raw EBSD scans, with six-row sample and 42-row sample--section MRD summaries for later manuscript-section selection.

**Architecture:** A MATLAB entry point reuses the registered catalog and loader, estimates one pixel-weighted ODF per scan with `6/mmm` and `SS = 1`, evaluates `phi2 = 0:10:60 deg`, writes both summaries, derives one global MRD limit from all 42 section maxima, then renders the matrix. It does not select manuscript representative sections.

**Tech Stack:** MATLAB R2025a, MTEX 6.1.1, project EBSD helpers, MATLAB `exportgraphics`.

## Global Constraints

- Use only six raw CTF files from `comprehensive_ebsd_catalog` and indexed `Ti-Hex` pixels with pixel weighting.
- Use `6/mmm` crystal symmetry, the `622` rotational fundamental zone, and bar-specimen `SS = 1`; never apply `222`.
- Use one 5 deg De la Vallee Poussin kernel, 5 deg evaluation resolution, and `phi2 = 0:10:60 deg` for every sample.
- Use one MRD range from zero to the maximum across all 42 section maxima and a perceptually ordered non-rainbow map.
- Do not modify raw or denoised CTF files.
- Export `odf_diameter_full_sections.png`, `odf_diameter_full_sections.pdf`, `odf_diameter_summary.csv` (six rows), and `odf_diameter_section_summary.csv` (42 rows) under `results/mtex_odf_diameter_montage`.
- Select manuscript angles only after checking the full diagnostic and the 42-row CSV. A later caption must use `selected phi2 sections` and list their angles.

---

### Task 1: Define the output contract with a failing integration test

**Files:**

- Create: `tools/mtex/test_generate_odf_diameter_montage.m`
- Test: `tools/mtex/test_generate_odf_diameter_montage.m`

**Interface:** `[sampleSummary, sectionSummary] = generate_odf_diameter_montage(scanRoot, outputRoot)`.

- [ ] **Step 1: Write the failing test**

```matlab
function test_generate_odf_diameter_montage(scanRoot, outputRoot)
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded for this test.");
if isfolder(outputRoot), rmdir(outputRoot,"s"); end
[sampleSummary, sectionSummary] = generate_odf_diameter_montage(scanRoot, outputRoot);
assert(height(sampleSummary) == 6);
assert(height(sectionSummary) == 42);
assert(isequal(sampleSummary.sample,["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"]));
assert(isequal(sampleSummary.diameter_mm,[7;6.48;6.02;5.6;5.25;5]));
assert(isequal(sectionSummary.phi2_deg(1:7),(0:10:60)'));
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
assert(height(readtable(fullfile(outputRoot,"odf_diameter_summary.csv"))) == 6);
assert(height(readtable(fullfile(outputRoot,"odf_diameter_section_summary.csv"))) == 42);
end
```

- [ ] **Step 2: Run RED**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "startup_mtex('noMenu'); addpath('tools/mtex'); test_generate_odf_diameter_montage(string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')), string(fullfile(pwd,'.codex_tmp','odf-montage-red')));"
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
- Test: `tools/mtex/test_generate_odf_diameter_montage.m`

**Interface:** Returns `sampleSummary` (six rows) and `sectionSummary` (42 rows), and writes their CSV counterparts.

- [ ] **Step 1: Implement configuration and catalog selection**

Set `diameterOrder = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"]`, `phi2Deg = (0:10:60)'`, `kernelHalfwidthDeg = 5`, and `gridResolutionDeg = 5`. Instantiate `alphaTiCS = crystalSymmetry("6/mmm")` and `barSS = specimenSymmetry("1")`; assert those values and the `622` rotational fundamental zone. Select only raw rows from `comprehensive_ebsd_catalog(scanRoot)`, assert six rows in `diameterOrder`, and reject `specimenSymmetry("222")` in every helper.

- [ ] **Step 2: Implement one ODF and seven maxima per sample**

For each catalog row, load `Ti-Hex`, calculate pixel-weighted density with `SO3DeLaValleePoussinKernel("halfwidth",5*degree)`, `CS = alphaTiCS`, and `SS = barSS`, then normalize by `double(mean(rbfOdf))`. Evaluate each `phi2Deg` in the `622` rotational fundamental zone at `5*degree`; store finite positive results in `sectionMaximumMrd(6,7)`. Set `globalMaximumMrd = max(sectionMaximumMrd,[],"all")`.

- [ ] **Step 3: Implement exact table schemas and CSV writes**

`sampleSummary` columns are `sample`, `diameter_mm`, `cold_reduction_percent`, `input_path`, `valid_ti_hex_orientation_count`, `crystal_symmetry`, `rotational_fundamental_zone`, `specimen_symmetry`, `kernel_halfwidth_deg`, `grid_resolution_deg`, `maximum_section_mrd`, and `global_color_limit_max_mrd`.

`sectionSummary` has those leading sample metadata columns plus `phi2_deg`, `section_maximum_mrd`, and `global_color_limit_max_mrd`. Build the 42 rows using `repelem(catalog.sample,7)`, `repmat(phi2Deg,6,1)`, and `reshape(sectionMaximumMrd.',[],1)`, so records are ordered by diameter then phi2. Write the tables with:

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

Add `imfinfo(fullfile(outputRoot,"odf_diameter_full_sections.png"))` checks for positive dimensions and `Height > Width/2`, then run the Task 1 command using `.codex_tmp/odf-montage-render-red`. Expected: failure because the renderer has not written the PNG/PDF.

- [ ] **Step 2: Implement the matrix renderer**

Create `tiledlayout(6,7)`, rows in registered diameter order and columns in `phi2Deg` order. Render one section per tile using:

```matlab
plotSection(odfModels{scanIndex},"phi2",phi2Deg(sectionIndex)*degree, ...
  "contourf","silent","resolution",gridResolutionDeg*degree, ...
  "colorRange",[0 globalMaximumMrd]);
mtexColorMap parula
```

Preserve the `622` and `SS = 1` assertions, add diameter/cold-reduction labels in the first column and phi2 headings in the first row, then attach one east-side `ODF intensity (MRD)` colorbar with `[0 globalMaximumMrd]`. Export `odf_diameter_full_sections.png` at 600 dpi and `odf_diameter_full_sections.pdf` with `exportgraphics`; assert both are non-empty.

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

Run `generate_odf_diameter_montage` in a fresh MATLAB process with the raw scans directory and `results/mtex_odf_diameter_montage`; expect four non-empty artifacts.

- [ ] **Step 2: Verify acceptance criteria**

Use `Import-Csv` to assert six sample rows, 42 section rows, unique values `6/mmm`, `622`, and `1` for the three symmetry columns, and exactly one global color-limit value equal to the maximum `section_maximum_mrd`.

- [ ] **Step 3: Inspect before selecting manuscript angles**

Open `odf_diameter_full_sections.png` and verify six ordered diameter rows, seven `phi2 = 0:10:60 deg` columns, readable labels, one shared colorbar, and no blank or clipped cells. Compare visible peaks with the 42-row summary before deciding any manuscript subset. The caption of any later manuscript figure must call the result `selected phi2 sections` and list the angles.

- [ ] **Step 4: Commit final artifacts**

```powershell
git add -- results/mtex_odf_diameter_montage
git commit -m "results: add multi-diameter ODF diagnostic"
```
