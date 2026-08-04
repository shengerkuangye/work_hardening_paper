# Six-State Alpha-Ti C-Axis Pole-Figure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate and verify one publication-ready 2x3 montage of pixel-weighted alpha-Ti `{0001}` pole figures for the six raw cold-reduction states.

**Architecture:** A focused MTEX generator loads the six raw scans through the registered catalog, calculates one normalized harmonic ODF per scan from all valid Ti-Hex points, determines one shared `{0001}` pole-density maximum, renders the montage, and writes auditable metadata. One integration test exercises the real scans and validates the scientific contract and exported files.

**Tech Stack:** MATLAB R2025a, MTEX 6.1.1, raw EDAX CTF scans, MATLAB `exportgraphics`.

## Global Constraints

- Use six raw scans from `comprehensive_ebsd_catalog`, ordered by cold reduction.
- Use valid `Ti-Hex` points with one equal weight per sampling point.
- Use a de la Vallee Poussin kernel with 5 degree halfwidth and a 5 degree PF evaluation grid.
- Plot antipodal `{0001}` equal-area upper-hemisphere pole figures with AD = X horizontal, TD/RD = Y, and ND = Z.
- Use one shared MRD range across all six panels and a perceptually ordered `parula` colormap.
- Do not modify raw CTF files or overwrite the six existing single-panel PF images.
- Run MATLAB hidden with `-nodesktop -nosplash -noFigureWindows -batch`; clean only a verified run-owned orphan `MATLABWindow.exe` after the main process exits.

---

### Task 1: Real-data montage contract and generator

**Files:**
- Create: `tools/mtex/test_generate_c_axis_pf_six_state_montage.m`
- Create: `tools/mtex/generate_c_axis_pf_six_state_montage.m`

**Interfaces:**
- Consumes: `generate_c_axis_pf_six_state_montage(scanRoot,outputDir)` with two scalar string paths.
- Produces: a six-row metadata table and the PNG, PDF, and CSV paths defined in the design specification.

- [ ] **Step 1: Write the failing integration test**

Create `tools/mtex/test_generate_c_axis_pf_six_state_montage.m` with a real-data test that calls the wished-for generator and asserts:

```matlab
function test_generate_c_axis_pf_six_state_montage(scanRoot)
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")),"MTEX must be loaded.");
outputDir = string(tempname);
mkdir(outputDir);
cleanupOutput = onCleanup(@() rmdir(outputDir,"s"));

metadata = generate_c_axis_pf_six_state_montage(scanRoot,outputDir);
assert(height(metadata) == 6);
assert(isequal(metadata.cold_reduction_percent, ...
  [0;14.31;26.04;36;43.75;48.98]));
assert(all(metadata.variant == "raw"));
assert(all(metadata.valid_ti_hex_point_count > 0));
assert(all(isfinite(metadata.panel_max_mrd) & metadata.panel_max_mrd > 0));
assert(all(metadata.shared_max_mrd == metadata.shared_max_mrd(1)));
assert(metadata.shared_max_mrd(1) >= max(metadata.panel_max_mrd));
assert(all(metadata.kernel_halfwidth_deg == 5));
assert(all(metadata.grid_resolution_deg == 5));
assert(all(metadata.weighting == "pixel_equal"));

pngPath = fullfile(outputDir,"c_axis_pf_six_state_montage.png");
pdfPath = fullfile(outputDir,"c_axis_pf_six_state_montage.pdf");
csvPath = fullfile(outputDir,"c_axis_pf_six_state_montage_metadata.csv");
assert(isfile(pngPath) && isfile(pdfPath) && isfile(csvPath));
info = imfinfo(pngPath);
assert(info.Width > info.Height && info.Width >= 6000);
[xDpi,yDpi] = image_resolution_dpi(info);
assert(abs(xDpi - 600) < 1 && abs(yDpi - 600) < 1);
pdfBytes = fileread(pdfPath);
assert(startsWith(pdfBytes,"%PDF-"));
assert(height(readtable(csvPath,"TextType","string")) == 6);

clear cleanupOutput
fprintf("test_generate_c_axis_pf_six_state_montage passed\n");
end

function [xDpi,yDpi] = image_resolution_dpi(info)
if strcmpi(info.ResolutionUnit,"meter")
  xDpi = info.XResolution * 0.0254;
  yDpi = info.YResolution * 0.0254;
else
  assert(strcmpi(info.ResolutionUnit,"inch"));
  xDpi = info.XResolution;
  yDpi = info.YResolution;
end
end
```

- [ ] **Step 2: Run the test and verify RED**

Run MATLAB hidden from the repository root with:

```matlab
startup_mtex("noMenu");
addpath(fullfile(pwd,"tools","mtex"));
scanRoot = string(fullfile(pwd,"data","ebsd_kpl_250221_7_df","scans"));
test_generate_c_axis_pf_six_state_montage(scanRoot);
```

Expected: FAIL because `generate_c_axis_pf_six_state_montage` is undefined.

- [ ] **Step 3: Implement the minimal generator**

Create `tools/mtex/generate_c_axis_pf_six_state_montage.m` with these exact responsibilities:

```matlab
function metadata = generate_c_axis_pf_six_state_montage(scanRoot,outputDir)
arguments
  scanRoot (1,1) string
  outputDir (1,1) string
end
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")),"MTEX must be loaded.");
if ~isfolder(outputDir), mkdir(outputDir); end

catalog = comprehensive_ebsd_catalog(scanRoot);
catalog = catalog(catalog.variant == "raw",:);
[~,order] = sort(catalog.cold_reduction_percent);
catalog = catalog(order,:);
assert(height(catalog) == 6);

kernelHalfwidthDeg = 5;
gridResolutionDeg = 5;
kernel = SO3DeLaValleePoussinKernel( ...
  "halfwidth",kernelHalfwidthDeg*degree);
odfs = cell(6,1);
pointCounts = zeros(6,1);
panelMaximum = zeros(6,1);
inputPaths = strings(6,1);

for scanIndex = 1:6
  [ebsdFull,~] = load_comprehensive_ebsd_scan(catalog(scanIndex,:));
  tiEbsd = ebsdFull("Ti-Hex");
  assert(~isempty(tiEbsd));
  pointCounts(scanIndex) = length(tiEbsd);
  inputPaths(scanIndex) = replace(fullfile( ...
    catalog.folder(scanIndex),catalog.input_file(scanIndex)),filesep,"/");
  rbfOdf = calcDensity(tiEbsd.orientations,"kernel",kernel, ...
    "weights",ones(length(tiEbsd),1),"silent");
  rbfOdf = rbfOdf / double(mean(rbfOdf));
  odfs{scanIndex} = SO3FunHarmonic(rbfOdf,"bandwidth",kernel.bandwidth);
  cAxis = Miller(0,0,0,1,tiEbsd.CS);
  density = calcPDF(odfs{scanIndex},cAxis,[],"antipodal");
  panelMaximum(scanIndex) = double(max(density, ...
    "resolution",gridResolutionDeg*degree));
end
sharedMaximum = max(panelMaximum);

render_montage(odfs,catalog,sharedMaximum,gridResolutionDeg, ...
  fullfile(outputDir,"c_axis_pf_six_state_montage.png"), ...
  fullfile(outputDir,"c_axis_pf_six_state_montage.pdf"));

sample = catalog.sample;
diameter_mm = catalog.diameter_mm;
cold_reduction_percent = catalog.cold_reduction_percent;
variant = catalog.variant;
input_path = inputPaths;
valid_ti_hex_point_count = pointCounts;
panel_max_mrd = panelMaximum;
shared_max_mrd = repmat(sharedMaximum,6,1);
kernel_halfwidth_deg = repmat(kernelHalfwidthDeg,6,1);
grid_resolution_deg = repmat(gridResolutionDeg,6,1);
pole = repmat("{0001}",6,1);
projection = repmat("equal-area upper hemisphere; antipodal",6,1);
specimen_axes = repmat("X=AD; Y=TD/RD; Z=ND",6,1);
weighting = repmat("pixel_equal",6,1);
metadata = table(sample,diameter_mm,cold_reduction_percent,variant, ...
  input_path,valid_ti_hex_point_count,panel_max_mrd,shared_max_mrd, ...
  kernel_halfwidth_deg,grid_resolution_deg,pole,projection, ...
  specimen_axes,weighting);
writetable(metadata,fullfile(outputDir, ...
  "c_axis_pf_six_state_montage_metadata.csv"));
end
```

Implement `render_montage` as a local function that:

```matlab
function render_montage(odfs,catalog,colorMaximum,gridResolutionDeg, ...
    pngPath,pdfPath)
previousAnnotations = getMTEXpref("pfAnnotations");
restoreAnnotations = onCleanup(@() setMTEXpref( ...
  "pfAnnotations",previousAnnotations));
setMTEXpref("pfAnnotations",@(varargin) text( ...
  [vector3d.X,vector3d.Y,vector3d.Z],{"AD","TD/RD","ND"}, ...
  "BackgroundColor","white","tag","axesLabels",varargin{:}));

figureHandle = figure("Visible","off","Color","white", ...
  "Units","inches","Position",[0.25 0.25 12 7.6]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle,2,3,"Padding","compact", ...
  "TileSpacing","compact");
colormap(figureHandle,parula(256));
for panelIndex = 1:6
  axesHandle = nexttile(layout);
  cAxis = Miller(0,0,0,1,odfs{panelIndex}.CS);
  plotPDF(odfs{panelIndex},cAxis,"antipodal","earea", ...
    "contourf","silent","parent",axesHandle, ...
    "resolution",gridResolutionDeg*degree, ...
    "colorRange",[0 colorMaximum]);
  title(axesHandle,sprintf("%.2f mm | %.2f%%", ...
    catalog.diameter_mm(panelIndex), ...
    catalog.cold_reduction_percent(panelIndex)), ...
    "Interpreter","none","FontWeight","normal");
end
colorbarHandle = colorbar(axesHandle,"eastoutside");
colorbarHandle.Layout.Tile = "east";
colorbarHandle.Label.String = "MRD";
title(layout,"Pixel-weighted alpha-Ti {0001} c-axis pole figures", ...
  "FontWeight","bold");
exportgraphics(figureHandle,pngPath,"Resolution",600, ...
  "BackgroundColor","white");
exportgraphics(figureHandle,pdfPath,"ContentType","vector", ...
  "BackgroundColor","white");
clear cleanupFigure restoreAnnotations
end
```

- [ ] **Step 4: Run the test and verify GREEN**

Run the same hidden MATLAB command from Step 2.

Expected: `test_generate_c_axis_pf_six_state_montage passed` and exit code 0.

- [ ] **Step 5: Run MATLAB Code Analyzer**

Run:

```matlab
files = { ...
  "tools/mtex/generate_c_axis_pf_six_state_montage.m", ...
  "tools/mtex/test_generate_c_axis_pf_six_state_montage.m"};
for k = 1:numel(files)
  issues = checkcode(files{k},"-id");
  assert(isempty(issues),files{k});
end
disp("C_AXIS_PF_CODECHECK_OK");
```

Expected: `C_AXIS_PF_CODECHECK_OK` and exit code 0.

- [ ] **Step 6: Commit the tested generator**

```powershell
git add -- tools/mtex/generate_c_axis_pf_six_state_montage.m tools/mtex/test_generate_c_axis_pf_six_state_montage.m
git commit -m "feat: add six-state c-axis pole figure"
```

---

### Task 2: Production export and visual verification

**Files:**
- Create: `results/mtex_c_axis_pf/c_axis_pf_six_state_montage.png`
- Create: `results/mtex_c_axis_pf/c_axis_pf_six_state_montage.pdf`
- Create: `results/mtex_c_axis_pf/c_axis_pf_six_state_montage_metadata.csv`

**Interfaces:**
- Consumes: the tested `generate_c_axis_pf_six_state_montage` function.
- Produces: the three final figure artifacts for manuscript use and audit.

- [ ] **Step 1: Generate production artifacts**

Run MATLAB hidden with:

```matlab
startup_mtex("noMenu");
addpath(fullfile(pwd,"tools","mtex"));
scanRoot = string(fullfile(pwd,"data","ebsd_kpl_250221_7_df","scans"));
outputDir = string(fullfile(pwd,"results","mtex_c_axis_pf"));
metadata = generate_c_axis_pf_six_state_montage(scanRoot,outputDir);
assert(height(metadata) == 6);
disp("C_AXIS_PF_PRODUCTION_OK");
```

Expected: `C_AXIS_PF_PRODUCTION_OK` and exit code 0.

- [ ] **Step 2: Verify output metadata and file properties**

Run:

```powershell
$out = 'results/mtex_c_axis_pf'
$png = Join-Path $out 'c_axis_pf_six_state_montage.png'
$pdf = Join-Path $out 'c_axis_pf_six_state_montage.pdf'
$csv = Join-Path $out 'c_axis_pf_six_state_montage_metadata.csv'
$rows = Import-Csv -Encoding UTF8 $csv
if ($rows.Count -ne 6) { throw 'Expected six metadata rows.' }
$actual = @($rows | ForEach-Object {[double]$_.cold_reduction_percent})
$expected = @(0,14.31,26.04,36,43.75,48.98)
for ($i=0; $i -lt 6; $i++) {
  if ([math]::Abs($actual[$i]-$expected[$i]) -gt 1e-8) {
    throw "Reduction order mismatch at row $i."
  }
}
Add-Type -AssemblyName System.Drawing
$image = [System.Drawing.Image]::FromFile((Resolve-Path $png))
try {
  if ($image.Width -le $image.Height -or $image.Width -lt 6000) {
    throw 'PNG dimensions do not satisfy the landscape contract.'
  }
  if ([math]::Abs($image.HorizontalResolution-600) -ge 1 -or
      [math]::Abs($image.VerticalResolution-600) -ge 1) {
    throw 'PNG resolution is not approximately 600 dpi.'
  }
} finally {
  $image.Dispose()
}
if ((Get-Item $pdf).Length -le 0) { throw 'PDF is empty.' }
$header = [System.Text.Encoding]::ASCII.GetString(
  [System.IO.File]::ReadAllBytes((Resolve-Path $pdf))[0..4])
if ($header -ne '%PDF-') { throw 'Invalid PDF header.' }
'PF_ARTIFACT_VERIFICATION_OK'
```

Expected: `PF_ARTIFACT_VERIFICATION_OK`.

- [ ] **Step 3: Visually inspect the PNG**

Open the full-resolution PNG and confirm:

- all six pole circles are complete;
- panel order and labels are correct;
- AD is horizontal in every panel;
- AD, TD/RD, and ND labels are readable;
- only one shared MRD colorbar appears;
- titles, contours, and colorbar are not clipped;
- all panels use the same color range.

- [ ] **Step 4: Run final focused verification**

Rerun the integration test, Code Analyzer, artifact assertions, `git diff --check`, and process inventory in one final verification pass.

Expected: no test failures, no Code Analyzer issues, no whitespace errors, and no run-owned MATLAB process left behind.

- [ ] **Step 5: Commit production artifacts**

```powershell
git add -- results/mtex_c_axis_pf/c_axis_pf_six_state_montage.png results/mtex_c_axis_pf/c_axis_pf_six_state_montage.pdf results/mtex_c_axis_pf/c_axis_pf_six_state_montage_metadata.csv docs/superpowers/plans/2026-08-04-six-state-c-axis-pole-figure.md
git commit -m "fig: add six-state c-axis pole figure montage"
```
