# Selected-phi2 ODF Figure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate a 6-row by 3-column alpha-Ti ODF figure at `phi2 = 0, 30, and 60 deg`, annotate each row with its formal MTEX ODF maximum, and write reproducible global-maximum and selected-section peak-position tables.

**Architecture:** A calculation function reuses the registered EBSD catalog, loader, and guarded normalization helper to return six harmonic ODFs, formal MTEX maxima, and selected-section peak tables without writing files. A focused generator calls that tested calculation interface, writes two CSVs, clips nonphysical negative values only on the contour display grid, and renders six identical three-section row rasters into one shared-scale figure.

**Tech Stack:** MATLAB R2025a, MTEX 6.1.1, existing project EBSD helpers, MATLAB `exportgraphics`.

## Global Constraints

- Preserve the existing seven-section diagnostic and its four outputs.
- Use only the six registered raw CTF scans and indexed `Ti-Hex` pixels with pixel weighting.
- Preserve imported `6/mmm`, proper rotational group `622`, and bar-specimen `SS = 1`; reject `222`.
- Use a 5 deg De la Vallee Poussin kernel for every sample.
- Use MTEX `max(odf,"resolution",1*degree)` for the formal ODF maximum.
- Use selected sections exactly `[0,30,60] deg`.
- Evaluate selected-section peak positions on `phi1 = 0:1:359 deg` and `Phi = 0:1:90 deg`.
- Use one color range from zero to the largest formal ODF maximum among all six samples.
- Clip negative contour-grid values to `0 MRD` for rendering only; do not change ODFs or CSV values.
- Export PNG and image-content PDF at approximately 600 dpi.
- Call the result **selected `phi2` sections**, not a complete ODF.

---

### Task 1: Define the selected-section output contract

**Files:**

- Create: `tools/mtex/test_generate_odf_selected_phi2_figure.m`

**Interfaces:**

- Consumes: `scanRoot (1,1) string`.
- Produces: an executable integration-test entry point
  `test_generate_odf_selected_phi2_figure(scanRoot)`.
- Defines the generator interface
  `[sampleSummary,peakSummary] = generate_odf_selected_phi2_figure(scanRoot,outputRoot)`.

- [ ] **Step 1: Write the failing integration test**

Create `tools/mtex/test_generate_odf_selected_phi2_figure.m`:

```matlab
function test_generate_odf_selected_phi2_figure(scanRoot)
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded for this test.");

outputRoot = string(tempname);
mkdir(outputRoot);
cleanupOutput = onCleanup(@() remove_test_output(outputRoot));
[sampleSummary,peakSummary] = ...
  generate_odf_selected_phi2_figure(scanRoot,outputRoot);

sampleColumns = ["sample","diameter_mm","cold_reduction_percent", ...
  "input_path","valid_ti_hex_orientation_count","crystal_symmetry", ...
  "rotational_fundamental_zone","specimen_symmetry", ...
  "kernel_halfwidth_deg","odf_maximum_resolution_deg", ...
  "odf_maximum_mrd","phi1_max_deg","Phi_max_deg","phi2_max_deg", ...
  "global_color_limit_max_mrd"];
peakColumns = ["sample","diameter_mm","cold_reduction_percent", ...
  "phi2_deg","section_peak_mrd","phi1_peak_deg","Phi_peak_deg", ...
  "section_grid_resolution_deg","global_color_limit_max_mrd"];

assert(isequal(string(sampleSummary.Properties.VariableNames),sampleColumns));
assert(isequal(string(peakSummary.Properties.VariableNames),peakColumns));
assert(height(sampleSummary) == 6);
assert(height(peakSummary) == 18);
assert(isequal(sampleSummary.sample, ...
  ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"]));
assert(isequal(sampleSummary.cold_reduction_percent, ...
  [0;14.31;26.04;36;43.75;48.98]));
assert(all(sampleSummary.crystal_symmetry == "6/mmm"));
assert(all(sampleSummary.rotational_fundamental_zone == "622"));
assert(all(sampleSummary.specimen_symmetry == "1"));
assert(all(sampleSummary.kernel_halfwidth_deg == 5));
assert(all(sampleSummary.odf_maximum_resolution_deg == 1));
assert(all(isfinite(sampleSummary.odf_maximum_mrd) & ...
  sampleSummary.odf_maximum_mrd > 0));
assert(all(sampleSummary.phi1_max_deg >= 0 & ...
  sampleSummary.phi1_max_deg <= 360));
assert(all(sampleSummary.Phi_max_deg >= 0 & ...
  sampleSummary.Phi_max_deg <= 180));
assert(all(sampleSummary.phi2_max_deg >= 0 & ...
  sampleSummary.phi2_max_deg <= 360));

for sampleIndex = 1:6
  rows = (sampleIndex - 1) * 3 + (1:3);
  assert(isequal(peakSummary.sample(rows), ...
    repmat(sampleSummary.sample(sampleIndex),3,1)));
  assert(isequal(peakSummary.phi2_deg(rows),[0;30;60]));
end
assert(all(isfinite(peakSummary.section_peak_mrd) & ...
  peakSummary.section_peak_mrd > 0));
assert(all(peakSummary.phi1_peak_deg >= 0 & ...
  peakSummary.phi1_peak_deg < 360));
assert(all(peakSummary.Phi_peak_deg >= 0 & ...
  peakSummary.Phi_peak_deg <= 90));
assert(all(peakSummary.section_grid_resolution_deg == 1));
globalLimit = max(sampleSummary.odf_maximum_mrd);
assert(all(sampleSummary.global_color_limit_max_mrd == globalLimit));
assert(all(peakSummary.global_color_limit_max_mrd == globalLimit));

sampleCsv = fullfile(outputRoot,"odf_selected_phi2_summary.csv");
peakCsv = fullfile(outputRoot,"odf_selected_phi2_peak_positions.csv");
pngPath = fullfile(outputRoot,"odf_selected_phi2_sections.png");
pdfPath = fullfile(outputRoot,"odf_selected_phi2_sections.pdf");
assert(isfile(sampleCsv) && isfile(peakCsv));
assert(isfile(pngPath) && isfile(pdfPath));
assert(height(readtable(sampleCsv,"TextType","string")) == 6);
assert(height(readtable(peakCsv,"TextType","string")) == 18);

imageInfo = imfinfo(pngPath);
assert(imageInfo.Width > imageInfo.Height);
assert(imageInfo.Width >= 2400 && imageInfo.Height >= 1800);
[xDpi,yDpi] = image_resolution_dpi(imageInfo);
assert(abs(xDpi - 600) < 1 && abs(yDpi - 600) < 1);
imageData = imread(pngPath);
assert(mean(any(imageData < 250,3),"all") >= 0.18);

pdfInfo = dir(pdfPath);
assert(pdfInfo.bytes > 0);
pdfText = fileread(pdfPath);
assert(startsWith(pdfText,"%PDF-"));
imageTokens = regexp(pdfText, ...
  '/Subtype /Image\s+/Width ([0-9]+)\s+/Height ([0-9]+)',"tokens");
mediaBox = regexp(pdfText, ...
  '/MediaBox \[0 0 ([0-9.]+) ([0-9.]+)\]',"tokens","once");
assert(~isempty(imageTokens) && ~isempty(mediaBox));
pixelArea = sum(cellfun(@(token) ...
  str2double(token{1}) * str2double(token{2}),imageTokens));
pageArea = (str2double(mediaBox{1}) / 72) * ...
  (str2double(mediaBox{2}) / 72);
assert(sqrt(pixelArea / pageArea) >= 590);
clear cleanupOutput
end

function [xDpi,yDpi] = image_resolution_dpi(imageInfo)
if strcmpi(imageInfo.ResolutionUnit,"meter")
  xDpi = imageInfo.XResolution * 0.0254;
  yDpi = imageInfo.YResolution * 0.0254;
else
  assert(strcmpi(imageInfo.ResolutionUnit,"Inch"));
  xDpi = imageInfo.XResolution;
  yDpi = imageInfo.YResolution;
end
end

function remove_test_output(path)
if isfolder(path)
  rmdir(path,"s");
end
end
```

- [ ] **Step 2: Run RED**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "startup_mtex('noMenu'); addpath('tools/mtex'); test_generate_odf_selected_phi2_figure(string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')));"
```

Expected: MATLAB exits non-zero because
`generate_odf_selected_phi2_figure` is undefined.

- [ ] **Step 3: Commit the failing contract**

```powershell
git add -- tools/mtex/test_generate_odf_selected_phi2_figure.m
git commit -m "test: define selected-phi2 ODF figure contract"
```

### Task 2: Calculate formal ODF maxima and selected-section peaks

**Files:**

- Create: `tools/mtex/calculate_selected_phi2_odf_data.m`
- Create: `tools/mtex/test_calculate_selected_phi2_odf_data.m`

**Interfaces:**

- Consumes: `scanRoot (1,1) string`.
- Produces:
  `[sampleSummary,peakSummary,odfs,catalog] = calculate_selected_phi2_odf_data(scanRoot)`.
- `sampleSummary`: six rows with formal MTEX ODF maxima and maximum
  orientations.
- `peakSummary`: 18 rows with selected-section peak values and positions.
- `odfs`: six normalized harmonic ODF objects used by Task 3.
- `catalog`: six registered raw catalog rows in decreasing-diameter order.

- [ ] **Step 1: Write the failing calculation test**

Create `tools/mtex/test_calculate_selected_phi2_odf_data.m`:

```matlab
function test_calculate_selected_phi2_odf_data(scanRoot)
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded for this test.");
[sampleSummary,peakSummary,odfs,catalog] = ...
  calculate_selected_phi2_odf_data(scanRoot);
assert(height(sampleSummary) == 6 && height(peakSummary) == 18);
assert(numel(odfs) == 6 && height(catalog) == 6);
assert(isequal(sampleSummary.sample, ...
  ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"]));
assert(isequal(catalog.sample,sampleSummary.sample));
assert(isequal(sampleSummary.cold_reduction_percent, ...
  [0;14.31;26.04;36;43.75;48.98]));
assert(all(sampleSummary.crystal_symmetry == "6/mmm"));
assert(all(sampleSummary.rotational_fundamental_zone == "622"));
assert(all(sampleSummary.specimen_symmetry == "1"));
assert(all(sampleSummary.odf_maximum_resolution_deg == 1));
assert(all(isfinite(sampleSummary.odf_maximum_mrd) & ...
  sampleSummary.odf_maximum_mrd > 0));
for sampleIndex = 1:6
  rows = (sampleIndex - 1) * 3 + (1:3);
  assert(isequal(peakSummary.phi2_deg(rows),[0;30;60]));
  assert(string(odfs{sampleIndex}.CS.pointGroup) == "6/mmm");
  assert(string(odfs{sampleIndex}.CS.properGroup.pointGroup) == "622");
  assert(string(odfs{sampleIndex}.SS.pointGroup) == "1");
end
assert(all(isfinite(peakSummary.section_peak_mrd) & ...
  peakSummary.section_peak_mrd > 0));
assert(all(peakSummary.section_grid_resolution_deg == 1));
globalLimit = max(sampleSummary.odf_maximum_mrd);
assert(all(sampleSummary.global_color_limit_max_mrd == globalLimit));
assert(all(peakSummary.global_color_limit_max_mrd == globalLimit));
end
```

- [ ] **Step 2: Run the calculation test to verify RED**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "startup_mtex('noMenu'); addpath('tools/mtex'); test_calculate_selected_phi2_odf_data(string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')));"
```

Expected: MATLAB exits non-zero because
`calculate_selected_phi2_odf_data` is undefined.

- [ ] **Step 3: Implement catalog selection and ODF construction**

Create `tools/mtex/calculate_selected_phi2_odf_data.m` with:

```matlab
function [sampleSummary,peakSummary,odfs,catalog] = ...
  calculate_selected_phi2_odf_data(scanRoot)
%CALCULATE_SELECTED_PHI2_ODF_DATA Calculate alpha-Ti ODF peak data.

arguments
  scanRoot (1,1) string
end
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded.");

diameterOrder = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"];
selectedPhi2Deg = [0;30;60];
kernelHalfwidthDeg = 5;
maximumResolutionDeg = 1;
sectionGridResolutionDeg = 1;

catalog = comprehensive_ebsd_catalog(scanRoot);
catalog = catalog(catalog.variant == "raw",:);
[isRegistered,catalogOrder] = ismember(diameterOrder,catalog.sample);
assert(all(isRegistered) && numel(unique(catalogOrder)) == 6);
catalog = catalog(catalogOrder,:);
assert(height(catalog) == 6 && isequal(catalog.sample,diameterOrder));

kernel = SO3DeLaValleePoussinKernel( ...
  "halfwidth",kernelHalfwidthDeg * degree);
sampleCount = height(catalog);
sectionCount = numel(selectedPhi2Deg);
odfs = cell(sampleCount,1);
validCounts = zeros(sampleCount,1);
odfMaximumMrd = zeros(sampleCount,1);
maximumEulerDeg = zeros(sampleCount,3);
sectionPeakMrd = zeros(sampleCount,sectionCount);
sectionPeakPhi1Deg = zeros(sampleCount,sectionCount);
sectionPeakPhiDeg = zeros(sampleCount,sectionCount);
phi1Deg = 0:sectionGridResolutionDeg:359;
PhiDeg = 0:sectionGridResolutionDeg:90;
[phi1Grid,PhiGrid] = meshgrid(phi1Deg,PhiDeg);

for scanIndex = 1:sampleCount
  [ebsdFull,~] = load_comprehensive_ebsd_scan(catalog(scanIndex,:));
  orientations = ebsdFull("Ti-Hex").orientations;
  assert(~isempty(orientations));
  assert(string(orientations.CS.pointGroup) == "6/mmm");
  assert(string(orientations.CS.properGroup.pointGroup) == "622");
  assert(string(orientations.SS.pointGroup) == "1");
  validCounts(scanIndex) = numel(orientations);

  rbfOdf = calcDensity(orientations,"kernel",kernel, ...
    "weights",ones(numel(orientations),1),"silent");
  rbfOdf = normalize_positive_mean_density(rbfOdf);
  odf = SO3FunHarmonic(rbfOdf,"bandwidth",kernel.bandwidth);
  odf = normalize_positive_mean_density(odf);
  assert(string(odf.CS.pointGroup) == "6/mmm");
  assert(string(odf.CS.properGroup.pointGroup) == "622");
  assert(string(odf.SS.pointGroup) == "1");
  odfs{scanIndex} = odf;

  [maximumValue,maximumOrientation] = max(odf, ...
    "resolution",maximumResolutionDeg * degree);
  odfMaximumMrd(scanIndex) = double(maximumValue);
  [phi1Max,PhiMax,phi2Max] = Euler(maximumOrientation);
  maximumEulerDeg(scanIndex,:) = [ ...
    mod(phi1Max / degree,360),PhiMax / degree, ...
    mod(phi2Max / degree,360)];

  for sectionIndex = 1:sectionCount
    sectionOrientations = orientation.byEuler( ...
      phi1Grid(:) * degree,PhiGrid(:) * degree, ...
      selectedPhi2Deg(sectionIndex) * degree, ...
      orientations.CS,orientations.SS);
    values = real(eval(odf,sectionOrientations));
    assert(all(isfinite(values)));
    [sectionPeakMrd(scanIndex,sectionIndex),peakIndex] = max(values);
    sectionPeakPhi1Deg(scanIndex,sectionIndex) = phi1Grid(peakIndex);
    sectionPeakPhiDeg(scanIndex,sectionIndex) = PhiGrid(peakIndex);
  end
  clear ebsdFull orientations rbfOdf odf
end
assert(all(isfinite(odfMaximumMrd) & odfMaximumMrd > 0));
assert(all(isfinite(sectionPeakMrd) & sectionPeakMrd > 0,"all"));
globalMaximumMrd = max(odfMaximumMrd);
```

- [ ] **Step 4: Build and return the exact tables**

Append:

```matlab
sample = catalog.sample;
diameter_mm = catalog.diameter_mm;
cold_reduction_percent = catalog.cold_reduction_percent;
input_path = catalog.input_path;
valid_ti_hex_orientation_count = validCounts;
crystal_symmetry = repmat("6/mmm",sampleCount,1);
rotational_fundamental_zone = repmat("622",sampleCount,1);
specimen_symmetry = repmat("1",sampleCount,1);
kernel_halfwidth_deg = repmat(kernelHalfwidthDeg,sampleCount,1);
odf_maximum_resolution_deg = repmat(maximumResolutionDeg,sampleCount,1);
odf_maximum_mrd = odfMaximumMrd;
phi1_max_deg = maximumEulerDeg(:,1);
Phi_max_deg = maximumEulerDeg(:,2);
phi2_max_deg = maximumEulerDeg(:,3);
global_color_limit_max_mrd = repmat(globalMaximumMrd,sampleCount,1);
sampleSummary = table(sample,diameter_mm,cold_reduction_percent, ...
  input_path,valid_ti_hex_orientation_count,crystal_symmetry, ...
  rotational_fundamental_zone,specimen_symmetry,kernel_halfwidth_deg, ...
  odf_maximum_resolution_deg,odf_maximum_mrd,phi1_max_deg, ...
  Phi_max_deg,phi2_max_deg,global_color_limit_max_mrd);

sample = repelem(catalog.sample,sectionCount);
diameter_mm = repelem(catalog.diameter_mm,sectionCount);
cold_reduction_percent = repelem( ...
  catalog.cold_reduction_percent,sectionCount);
phi2_deg = repmat(selectedPhi2Deg,sampleCount,1);
section_peak_mrd = reshape(sectionPeakMrd.',[],1);
phi1_peak_deg = reshape(sectionPeakPhi1Deg.',[],1);
Phi_peak_deg = reshape(sectionPeakPhiDeg.',[],1);
section_grid_resolution_deg = repmat( ...
  sectionGridResolutionDeg,sampleCount * sectionCount,1);
global_color_limit_max_mrd = repmat( ...
  globalMaximumMrd,sampleCount * sectionCount,1);
peakSummary = table(sample,diameter_mm,cold_reduction_percent, ...
  phi2_deg,section_peak_mrd,phi1_peak_deg,Phi_peak_deg, ...
  section_grid_resolution_deg,global_color_limit_max_mrd);

assert(height(sampleSummary) == 6 && height(peakSummary) == 18);
assert(isequal(reshape(peakSummary.phi2_deg,3,[]), ...
  repmat(selectedPhi2Deg,1,6)));
end
```

- [ ] **Step 5: Run the focused calculation test to verify GREEN**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "startup_mtex('noMenu'); addpath('tools/mtex'); test_calculate_selected_phi2_odf_data(string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')));"
```

Expected: MATLAB exits 0 and all ordering, symmetry, formal-maximum,
selected-section, and shared-color-limit assertions pass without creating
output files.

- [ ] **Step 6: Commit the calculation interface**

```powershell
git add -- tools/mtex/calculate_selected_phi2_odf_data.m tools/mtex/test_calculate_selected_phi2_odf_data.m
git commit -m "feat: calculate selected-phi2 ODF maxima"
```

### Task 3: Render the 6-row by 3-column publication figure

**Files:**

- Create: `tools/mtex/generate_odf_selected_phi2_figure.m`
- Modify: `tools/mtex/test_generate_odf_selected_phi2_figure.m`

**Interfaces:**

- Consumes: six harmonic ODFs, six catalog rows, `[0;30;60]`, six formal
  maxima, one global color limit, PNG path, and PDF path.
- Produces: one 600 dpi PNG and one 600 dpi image-content PDF.

- [ ] **Step 1: Add the generator orchestration**

Create `tools/mtex/generate_odf_selected_phi2_figure.m`:

```matlab
function [sampleSummary,peakSummary] = ...
  generate_odf_selected_phi2_figure(scanRoot,outputRoot)
%GENERATE_ODF_SELECTED_PHI2_FIGURE Render selected alpha-Ti ODF sections.

arguments
  scanRoot (1,1) string
  outputRoot (1,1) string
end
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded.");
if ~isfolder(outputRoot), mkdir(outputRoot); end
[sampleSummary,peakSummary,odfs,catalog] = ...
  calculate_selected_phi2_odf_data(scanRoot);
writetable(sampleSummary, ...
  fullfile(outputRoot,"odf_selected_phi2_summary.csv"));
writetable(peakSummary, ...
  fullfile(outputRoot,"odf_selected_phi2_peak_positions.csv"));
selectedPhi2Deg = [0;30;60];
pngPath = fullfile(outputRoot,"odf_selected_phi2_sections.png");
pdfPath = fullfile(outputRoot,"odf_selected_phi2_sections.pdf");
render_selected_phi2_matrix(odfs,catalog,selectedPhi2Deg, ...
  sampleSummary.odf_maximum_mrd, ...
  sampleSummary.global_color_limit_max_mrd(1),pngPath,pdfPath);
end
```

- [ ] **Step 2: Add the row-raster renderer**

Add this local function after the generator:

```matlab
function render_selected_phi2_matrix(odfs,catalog,selectedPhi2Deg, ...
  odfMaximumMrd,globalMaximumMrd,pngPath,pdfPath)
sampleCount = height(catalog);
sectionCount = numel(selectedPhi2Deg);
assert(sampleCount == 6 && sectionCount == 3);
tempRoot = string(tempname);
mkdir(tempRoot);
cleanupTemp = onCleanup(@() remove_temp_directory(tempRoot));
rowPaths = strings(sampleCount,1);
rowSizes = zeros(sampleCount,2);
renderedSectionCount = 0;

for scanIndex = 1:sampleCount
  rowPaths(scanIndex) = fullfile(tempRoot, ...
    sprintf("selected_phi2_row_%d.png",scanIndex));
  rowFigure = figure("Visible","off","Color","white", ...
    "Units","pixels","Position",[50 50 1100 520]);
  cleanupRow = onCleanup(@() close(rowFigure));
  plotSection(odfs{scanIndex},"phi2",selectedPhi2Deg * degree, ...
    "contourf","silent","layout",[1 sectionCount], ...
    "resolution",5 * degree,"colorRange",[0 globalMaximumMrd]);
  sectionAxes = findall(rowFigure,"Type","axes");
  sectionAxes = sectionAxes(arrayfun(@(axisHandle) ...
    isappdata(axisHandle,"sphericalPlot"),sectionAxes));
  assert(numel(sectionAxes) == sectionCount);
  assert(all(arrayfun(@(axisHandle) ...
    ~isempty(axisHandle.Children),sectionAxes)));
  clip_negative_contour_zdata(sectionAxes);
  for sectionIndex = 1:sectionCount
    sphericalPlotHandle = getappdata( ...
      sectionAxes(sectionIndex),"sphericalPlot");
    sphericalPlotHandle.TR.String = "";
  end
  mtexColorMap parula
  drawnow;
  exportgraphics(rowFigure,char(rowPaths(scanIndex)), ...
    "Resolution",600,"BackgroundColor","white");
  rowInfo = imfinfo(rowPaths(scanIndex));
  [rowDpiX,rowDpiY] = image_resolution_dpi(rowInfo);
  assert(abs(rowDpiX - 600) < 1 && abs(rowDpiY - 600) < 1);
  rowSizes(scanIndex,:) = [rowInfo.Width,rowInfo.Height];
  renderedSectionCount = renderedSectionCount + numel(sectionAxes);
  clear cleanupRow
end
assert(renderedSectionCount == 18);
assert(size(unique(rowSizes,"rows"),1) == 1);
crop_row_images_to_common_margins(rowPaths,36);

rowImages = cellfun(@imread,cellstr(rowPaths),"UniformOutput",false);
rowWidths = cellfun(@(imageData) size(imageData,2),rowImages);
rowHeights = cellfun(@(imageData) size(imageData,1),rowImages);
assert(numel(unique(rowWidths)) == 1 && numel(unique(rowHeights)) == 1);
targetWidth = rowWidths(1);
rowHeight = rowHeights(1);
outputDpi = 600;
labelWidthPx = max(520,round(0.20 * targetWidth));
rowGapPx = max(12,round(0.025 * rowHeight));
titleHeightPx = max(90,round(0.14 * rowHeight));
columnTitleHeightPx = max(80,round(0.12 * rowHeight));
topGapPx = max(16,round(0.025 * rowHeight));
bottomPadPx = max(24,round(0.035 * rowHeight));
colorbarGapPx = max(50,round(0.014 * targetWidth));
colorbarWidthPx = max(70,round(0.018 * targetWidth));
colorbarLabelWidthPx = max(180,round(0.055 * targetWidth));
rightPadPx = max(30,round(0.008 * targetWidth));
rowsHeightPx = sampleCount * rowHeight + ...
  (sampleCount - 1) * rowGapPx;
topPadPx = titleHeightPx + columnTitleHeightPx + topGapPx;
canvasWidthPx = labelWidthPx + targetWidth + colorbarGapPx + ...
  colorbarWidthPx + colorbarLabelWidthPx + rightPadPx;
canvasHeightPx = bottomPadPx + rowsHeightPx + topPadPx;

matrixFigure = figure("Visible","off","Color","white", ...
  "Units","inches","Position",[0.25 0.25 ...
  canvasWidthPx / outputDpi,canvasHeightPx / outputDpi]);
cleanupMatrix = onCleanup(@() close(matrixFigure));
for scanIndex = 1:sampleCount
  rowBottomPx = bottomPadPx + ...
    (sampleCount - scanIndex) * (rowHeight + rowGapPx);
  rowAxes = axes(matrixFigure,"Units","normalized","Position",[ ...
    labelWidthPx / canvasWidthPx,rowBottomPx / canvasHeightPx, ...
    targetWidth / canvasWidthPx,rowHeight / canvasHeightPx]);
  image(rowAxes,rowImages{scanIndex});
  axis(rowAxes,"off");
  set(rowAxes,"DataAspectRatio",[1 1 1], ...
    "PlotBoxAspectRatio",[targetWidth rowHeight 1], ...
    "PositionConstraint","innerposition");
  annotation(matrixFigure,"textbox",[ ...
    0,rowBottomPx / canvasHeightPx, ...
    0.96 * labelWidthPx / canvasWidthPx,rowHeight / canvasHeightPx], ...
    "String",sprintf("%.2f mm\n%.2f%% reduction\nMax ODF = %.2f MRD", ...
    catalog.diameter_mm(scanIndex), ...
    catalog.cold_reduction_percent(scanIndex), ...
    odfMaximumMrd(scanIndex)), ...
    "HorizontalAlignment","right","VerticalAlignment","middle", ...
    "Interpreter","none","FontSize",8,"Color","black", ...
    "LineStyle","none","Margin",0);
end

rowsTopPx = bottomPadPx + rowsHeightPx;
columnWidthPx = targetWidth / sectionCount;
for sectionIndex = 1:sectionCount
  columnLeftPx = labelWidthPx + ...
    (sectionIndex - 1) * columnWidthPx;
  annotation(matrixFigure,"textbox",[ ...
    columnLeftPx / canvasWidthPx, ...
    (rowsTopPx + topGapPx) / canvasHeightPx, ...
    columnWidthPx / canvasWidthPx, ...
    columnTitleHeightPx / canvasHeightPx], ...
    "String",sprintf("\\phi_2 = %d^\\circ", ...
    selectedPhi2Deg(sectionIndex)), ...
    "HorizontalAlignment","center","VerticalAlignment","middle", ...
    "Interpreter","tex","FontSize",9,"Color","black", ...
    "LineStyle","none","Margin",0);
end
annotation(matrixFigure,"textbox",[ ...
  labelWidthPx / canvasWidthPx, ...
  (rowsTopPx + topGapPx + columnTitleHeightPx) / canvasHeightPx, ...
  targetWidth / canvasWidthPx,titleHeightPx / canvasHeightPx], ...
  "String","Selected ODF sections (SS = 1)", ...
  "HorizontalAlignment","center","VerticalAlignment","middle", ...
  "Interpreter","none","FontSize",11,"Color","black", ...
  "LineStyle","none","Margin",0);

colormap(matrixFigure,parula);
colorbarLeftPx = labelWidthPx + targetWidth + colorbarGapPx;
scaleAxes = axes(matrixFigure,"Units","normalized","Position",[ ...
  colorbarLeftPx / canvasWidthPx,bottomPadPx / canvasHeightPx, ...
  1 / canvasWidthPx,rowsHeightPx / canvasHeightPx],"Visible","off");
clim(scaleAxes,[0 globalMaximumMrd]);
colorbarHandle = colorbar(scaleAxes,"eastoutside");
colorbarHandle.Units = "normalized";
colorbarHandle.Position = [ ...
  colorbarLeftPx / canvasWidthPx,bottomPadPx / canvasHeightPx, ...
  colorbarWidthPx / canvasWidthPx,rowsHeightPx / canvasHeightPx];
colorbarHandle.Label.String = "ODF intensity (MRD)";
colorbarHandle.Limits = [0 globalMaximumMrd];
colorbarHandle.FontSize = 8;
drawnow;

exportgraphics(matrixFigure,char(pngPath),"Resolution",600, ...
  "BackgroundColor","white");
exportgraphics(matrixFigure,char(pdfPath), ...
  "ContentType","image","Resolution",600, ...
  "BackgroundColor","white");
clear cleanupMatrix cleanupTemp
end
```

- [ ] **Step 3: Add local image and cleanup helpers**

Append:

```matlab
function [xDpi,yDpi] = image_resolution_dpi(imageInfo)
if strcmpi(imageInfo.ResolutionUnit,"meter")
  xDpi = imageInfo.XResolution * 0.0254;
  yDpi = imageInfo.YResolution * 0.0254;
else
  assert(strcmpi(imageInfo.ResolutionUnit,"Inch"));
  xDpi = imageInfo.XResolution;
  yDpi = imageInfo.YResolution;
end
end

function crop_row_images_to_common_margins(paths,padding)
rowImages = cellfun(@imread,cellstr(paths),"UniformOutput",false);
referenceSize = size(rowImages{1});
assert(all(cellfun(@(imageData) ...
  isequal(size(imageData),referenceSize),rowImages)));
unionNonwhite = false(referenceSize(1),referenceSize(2));
for rowIndex = 1:numel(rowImages)
  unionNonwhite = unionNonwhite | any(rowImages{rowIndex} < 250,3);
end
[contentRows,contentColumns] = find(unionNonwhite);
assert(~isempty(contentRows));
firstRow = max(1,min(contentRows) - padding);
lastRow = min(referenceSize(1),max(contentRows) + padding);
firstColumn = max(1,min(contentColumns) - padding);
lastColumn = min(referenceSize(2),max(contentColumns) + padding);
for rowIndex = 1:numel(rowImages)
  imageData = rowImages{rowIndex};
  imwrite(imageData(firstRow:lastRow,firstColumn:lastColumn,:), ...
    paths(rowIndex));
end
end

function remove_temp_directory(path)
if isfolder(path)
  rmdir(path,"s");
end
end
```

- [ ] **Step 4: Run GREEN**

Run the Task 1 command.

Expected: MATLAB exits 0; the temporary output directory is deleted by the
test; all table, symmetry, selected-section, raster, PDF, and resolution
assertions pass.

- [ ] **Step 5: Commit the generator and renderer**

```powershell
git add -- tools/mtex/generate_odf_selected_phi2_figure.m tools/mtex/test_generate_odf_selected_phi2_figure.m
git commit -m "feat: render selected-phi2 ODF figure"
```

### Task 4: Generate and verify formal outputs

**Files:**

- Create: `results/mtex_odf_selected_phi2/odf_selected_phi2_sections.png`
- Create: `results/mtex_odf_selected_phi2/odf_selected_phi2_sections.pdf`
- Create: `results/mtex_odf_selected_phi2/odf_selected_phi2_summary.csv`
- Create: `results/mtex_odf_selected_phi2/odf_selected_phi2_peak_positions.csv`

**Interfaces:**

- Consumes: committed generator and six registered raw scans.
- Produces: the four formal derived outputs.

- [ ] **Step 1: Generate formal results in a fresh MATLAB process**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "startup_mtex('noMenu'); addpath('tools/mtex'); generate_odf_selected_phi2_figure(string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')), string(fullfile(pwd,'results','mtex_odf_selected_phi2')));"
```

Expected: exit 0 and four non-empty files.

- [ ] **Step 2: Verify numerical and file contracts**

```powershell
$root = 'results\mtex_odf_selected_phi2'
$sample = Import-Csv (Join-Path $root 'odf_selected_phi2_summary.csv')
$peak = Import-Csv (Join-Path $root 'odf_selected_phi2_peak_positions.csv')
if ($sample.Count -ne 6) { throw 'Expected six ODF maximum rows.' }
if ($peak.Count -ne 18) { throw 'Expected 18 selected-section peak rows.' }
if (($sample.sample -join ',') -ne '7d,6.48d,6.02d,5.6d,5.25d,5d') {
  throw 'Unexpected sample order.'
}
foreach ($name in $sample.sample) {
  $angles = @($peak | Where-Object sample -eq $name |
    ForEach-Object { [double]$_.phi2_deg })
  if (($angles -join ',') -ne '0,30,60') {
    throw "Unexpected selected sections for $name."
  }
}
$global = @($sample.global_color_limit_max_mrd | Sort-Object -Unique)
$maximum = [double](($sample.odf_maximum_mrd |
  Measure-Object -Maximum).Maximum)
if ($global.Count -ne 1 -or
    [math]::Abs(([double]$global[0]) - $maximum) -gt 1e-12) {
  throw 'Global color limit does not equal the formal ODF maximum.'
}
Get-Item (Join-Path $root '*') | Select-Object Name,Length
```

Expected: no exception; six and 18 rows; one global limit equal to the formal
maximum; four non-empty files.

- [ ] **Step 3: Inspect the final figure**

Open `odf_selected_phi2_sections.png` and verify:

- six ordered rows and three columns;
- column headings `phi2 = 0, 30, 60 deg`;
- left labels contain diameter, registered reduction, and `Max ODF`;
- black contours remain visible;
- no white contour holes;
- one shared far-right MRD colorbar;
- no text or panel clipping.

- [ ] **Step 4: Run the fresh completion gate**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "startup_mtex('noMenu'); addpath('tools/mtex'); scanRoot=string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')); test_comprehensive_ebsd_contract; test_comprehensive_ebsd_helpers(scanRoot); test_calculate_selected_phi2_odf_data(scanRoot); test_generate_odf_selected_phi2_figure(scanRoot);"
git diff --check
git status --short
```

Expected: MATLAB exits 0, both baseline tests print `passed`, the selected
ODF test exits without assertion errors, `git diff --check` is silent, and
only intentional result files are uncommitted.

- [ ] **Step 5: Commit formal outputs**

```powershell
git add -- results/mtex_odf_selected_phi2/odf_selected_phi2_sections.png results/mtex_odf_selected_phi2/odf_selected_phi2_sections.pdf results/mtex_odf_selected_phi2/odf_selected_phi2_summary.csv results/mtex_odf_selected_phi2/odf_selected_phi2_peak_positions.csv
git commit -m "results: add selected-phi2 ODF figure"
```
