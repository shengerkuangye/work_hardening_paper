# OIM-style Square ODF Peak Matrix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate a 6-diameter by 7-section square ODF matrix at `phi2 = 0:10:60 deg`, mark each displayed-region peak, and export auditable peak-orientation tables without imposing specimen symmetry.

**Architecture:** A calculation-only MATLAB function loads the six registered raw CTF scans, builds normalized harmonic `SS = 1` ODFs, evaluates square and complete section grids, and returns tables plus the 42 square grids. A separate generator renders those grids with MATLAB axes, writes the two CSV records, and exports 600 dpi PNG/PDF files. Calculation and rendering have separate assertion-based integration tests so failures can be isolated.

**Tech Stack:** MATLAB R2025a, MTEX 6.1.1, project EBSD catalog/loading helpers, MATLAB tables and graphics, PowerShell hidden-process orchestration, PNG/PDF/CSV outputs.

## Global Constraints

- Preserve the six raw CTF inputs; create only new source, tests, plan, and files under `results/mtex_odf_square_peak_matrix/`.
- Use samples in this exact order: `7d`, `6.48d`, `6.02d`, `5.6d`, `5.25d`, `5d`.
- Use indexed raw `Ti-Hex` pixels, pixel weighting, imported `6/mmm` / `622`, and `SS = 1`; reject `222`/`mmm`.
- Use a 5 deg De la Vallee Poussin kernel and normalize the harmonic ODF to unit mean MRD.
- Evaluate `phi2 = 0:10:60 deg` on a 1 deg grid.
- The plotted domain is `phi1 = 0:90 deg`, `Phi = 0:90 deg`; the audit domain is `phi1 = 0:359 deg`, `Phi = 0:90 deg`.
- A plotted peak is always described as a maximum within the displayed square region.
- Use one shared color upper limit equal to the largest of the 42 square-domain peaks.
- Run every MATLAB batch with `Start-Process -WindowStyle Hidden` and arguments `-nodesktop -nosplash -noFigureWindows -batch`.
- Before each MATLAB batch, inventory existing `matlab.exe`, `MATLAB.exe`, and `MATLABWindow.exe` PIDs. After exit, inventory again and terminate only a new, exact run-owned `MATLABWindow.exe` PID after verifying no run-owned launcher/main MATLAB PID remains. Never terminate a pre-existing process.
- Stage and commit only task-owned paths. Do not alter the existing full-`phi1` or selected-`phi2` ODF generators or their results.

Use this PowerShell helper in the same shell for every MATLAB command below:

```powershell
function Invoke-ProjectMatlabBatch {
  param(
    [Parameter(Mandatory)][string]$Batch,
    [Parameter(Mandatory)][string]$LogStem
  )
  $matlabExe = 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe'
  if (-not (Test-Path -LiteralPath $matlabExe)) {
    throw "MATLAB executable not found: $matlabExe"
  }
  $logRoot = Join-Path (Get-Location) '.codex_tmp\matlab_logs'
  New-Item -ItemType Directory -Force -Path $logRoot | Out-Null
  $stdoutPath = Join-Path $logRoot ($LogStem + '.stdout.txt')
  $stderrPath = Join-Path $logRoot ($LogStem + '.stderr.txt')
  $before = @(Get-Process -ErrorAction SilentlyContinue |
    Where-Object { $_.ProcessName -in @('matlab','MATLAB','MATLABWindow') } |
    Select-Object Id,ProcessName,StartTime)
  $beforePids = @($before.Id)
  $process = Start-Process -FilePath $matlabExe `
    -ArgumentList @(
      '-nodesktop',
      '-nosplash',
      '-noFigureWindows',
      '-batch',
      $Batch
    ) `
    -WindowStyle Hidden `
    -RedirectStandardOutput $stdoutPath `
    -RedirectStandardError $stderrPath `
    -PassThru `
    -Wait
  $exitCode = $process.ExitCode
  $after = @(Get-Process -ErrorAction SilentlyContinue |
    Where-Object { $_.ProcessName -in @('matlab','MATLAB','MATLABWindow') } |
    Select-Object Id,ProcessName,StartTime)
  $newProcesses = @($after | Where-Object { $_.Id -notin $beforePids })
  $runOwnedMain = @($newProcesses |
    Where-Object { $_.ProcessName -in @('matlab','MATLAB') })
  if ($runOwnedMain.Count -eq 0) {
    $runOwnedWindows = @($newProcesses |
      Where-Object { $_.ProcessName -eq 'MATLABWindow' })
    foreach ($window in $runOwnedWindows) {
      Stop-Process -Id $window.Id -Force
    }
  } elseif (@($newProcesses |
      Where-Object { $_.ProcessName -eq 'MATLABWindow' }).Count -gt 0) {
    throw 'Run-owned MATLAB main process still exists; refusing orphan cleanup.'
  }
  $finalInventory = @(Get-Process -ErrorAction SilentlyContinue |
    Where-Object { $_.ProcessName -in @('matlab','MATLAB','MATLABWindow') } |
    Select-Object Id,ProcessName,StartTime)
  [pscustomobject]@{
    ExitCode = $exitCode
    Before = $before
    After = $after
    Final = $finalInventory
    Stdout = $stdoutPath
    Stderr = $stderrPath
  } | Format-List
  Get-Content -LiteralPath $stdoutPath
  if ((Test-Path -LiteralPath $stderrPath) -and
      (Get-Item -LiteralPath $stderrPath).Length -gt 0) {
    Get-Content -LiteralPath $stderrPath
  }
  if ($exitCode -ne 0) {
    throw "MATLAB batch failed with exit code $exitCode"
  }
}
```

---

### Task 1: Calculate square and complete-section ODF peaks

**Files:**
- Create: `tools/mtex/test_calculate_square_odf_peak_data.m`
- Create: `tools/mtex/calculate_square_odf_peak_data.m`

**Interfaces:**
- Consumes: `comprehensive_ebsd_catalog(scanRoot)`, `load_comprehensive_ebsd_scan(catalogRow)`, and `normalize_positive_mean_density(odf)`.
- Produces: `[peakSummary,parameterSummary,odfs,catalog,squareSectionValues] = calculate_square_odf_peak_data(scanRoot)`.
- `peakSummary` is a 42-row table; `parameterSummary` is a one-row table; `odfs` is a 6-by-1 cell; `catalog` is the ordered 6-row raw catalog; `squareSectionValues` is a 6-by-7 cell of 91-by-91 MRD grids.

- [ ] **Step 1: Write the failing calculation contract test**

Create `tools/mtex/test_calculate_square_odf_peak_data.m`:

```matlab
function test_calculate_square_odf_peak_data(scanRoot)
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded for this test.");
originalConvention = getMTEXpref("EulerAngleConvention");
originalDirectory = pwd;
temporaryDirectory = string(tempname);
mkdir(temporaryDirectory);
cleanupState = onCleanup(@() restore_state( ...
  originalConvention,originalDirectory,temporaryDirectory));
setMTEXpref("EulerAngleConvention","Roe");
cd(temporaryDirectory);
[peaks,parameters,odfs,catalog,squareValues] = ...
  calculate_square_odf_peak_data(scanRoot);
temporaryContents = dir(temporaryDirectory);
temporaryContents = temporaryContents(~ismember( ...
  string({temporaryContents.name}),[".",".."]));
assert(isempty(temporaryContents), ...
  "Calculation function must not write files.");
cd(originalDirectory);

expectedSamples = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"];
expectedPhi2 = (0:10:60)';
assert(height(peaks) == 42 && height(parameters) == 1);
assert(numel(odfs) == 6 && height(catalog) == 6);
assert(isequal(catalog.sample,expectedSamples));
assert(all(catalog.variant == "raw"));
assert(isequal(size(squareValues),[6 7]));
assert(all(cellfun(@(values) isequal(size(values),[91 91]), ...
  squareValues),"all"));
assert(isequal(peaks.peak_id,compose("P%02d",(1:42)')));

for sampleIndex = 1:6
  rows = (sampleIndex - 1) * 7 + (1:7);
  assert(isequal(peaks.sample(rows), ...
    repmat(expectedSamples(sampleIndex),7,1)));
  assert(isequal(peaks.phi2_deg(rows),expectedPhi2));
  assert(string(odfs{sampleIndex}.CS.pointGroup) == "6/mmm");
  assert(string(odfs{sampleIndex}.CS.properGroup.pointGroup) == "622");
  assert(string(odfs{sampleIndex}.SS.pointGroup) == "1");
end

assert(all(isfinite(peaks.display_peak_mrd) & peaks.display_peak_mrd > 0));
assert(all(peaks.display_phi1_deg >= 0 & peaks.display_phi1_deg <= 90));
assert(all(peaks.display_Phi_deg >= 0 & peaks.display_Phi_deg <= 90));
assert(all(peaks.display_phi2_deg == peaks.phi2_deg));
assert(all(isfinite(peaks.complete_peak_mrd) & ...
  peaks.complete_peak_mrd > 0));
assert(all(peaks.complete_phi1_deg >= 0 & ...
  peaks.complete_phi1_deg <= 359));
assert(all(peaks.complete_Phi_deg >= 0 & ...
  peaks.complete_Phi_deg <= 90));
assert(all(peaks.complete_phi2_deg == peaks.phi2_deg));
assert(all(peaks.complete_peak_mrd + 1e-10 >= peaks.display_peak_mrd));
assert(all(abs(peaks.complete_minus_display_mrd - ...
  (peaks.complete_peak_mrd - peaks.display_peak_mrd)) <= 1e-10));
assert(isequal(peaks.complete_peak_inside_display, ...
  peaks.complete_phi1_deg <= 90));

axisAngles = [peaks.c_axis_to_ad_deg, ...
  peaks.c_axis_to_td_rd_deg,peaks.c_axis_to_nd_deg];
assert(all(isfinite(axisAngles),"all"));
assert(all(axisAngles >= 0 & axisAngles <= 90,"all"));
axisNames = ["AD","TD/RD","ND"];
for rowIndex = 1:height(peaks)
  [minimumAngle,minimumIndex] = min(axisAngles(rowIndex,:));
  expectedNote = sprintf("c-axis closest to %s (%.1f deg)", ...
    axisNames(minimumIndex),minimumAngle);
  assert(peaks.orientation_note(rowIndex) == expectedNote);
end

globalLimit = max(peaks.display_peak_mrd);
assert(all(peaks.global_color_limit_max_mrd == globalLimit));
assert(parameters.data_variant == "raw");
assert(parameters.phase == "Ti-Hex");
assert(parameters.weighting == "pixel");
assert(parameters.crystal_symmetry == "6/mmm");
assert(parameters.rotational_group == "622");
assert(parameters.specimen_symmetry == "1");
assert(parameters.kernel_halfwidth_deg == 5);
assert(parameters.grid_resolution_deg == 1);
assert(parameters.phi2_sections_deg == "0,10,20,30,40,50,60");
assert(parameters.display_phi1_range_deg == "0:1:90");
assert(parameters.display_Phi_range_deg == "0:1:90");
assert(parameters.complete_phi1_range_deg == "0:1:359");
assert(parameters.complete_Phi_range_deg == "0:1:90");
assert(getMTEXpref("EulerAngleConvention") == "Roe");
clear cleanupState
fprintf("test_calculate_square_odf_peak_data passed\n");
end

function restore_state(originalConvention,originalDirectory, ...
    temporaryDirectory)
setMTEXpref("EulerAngleConvention",originalConvention);
cd(originalDirectory);
if isfolder(temporaryDirectory)
  rmdir(temporaryDirectory,"s");
end
end
```

- [ ] **Step 2: Run the test and observe the expected failure**

```powershell
$batch = "addpath(fullfile(pwd,'tools','mtex')); startup_project_mtex; scanRoot=string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')); test_calculate_square_odf_peak_data(scanRoot);"
Invoke-ProjectMatlabBatch -Batch $batch -LogStem 'square_odf_calculation_red'
```

Expected: MATLAB exits nonzero because
`calculate_square_odf_peak_data` is undefined.

- [ ] **Step 3: Implement the calculation function**

Create `tools/mtex/calculate_square_odf_peak_data.m`:

```matlab
function [peakSummary,parameterSummary,odfs,catalog, ...
  squareSectionValues] = calculate_square_odf_peak_data(scanRoot)
%CALCULATE_SQUARE_ODF_PEAK_DATA Evaluate OIM-style square ODF sections.

arguments
  scanRoot (1,1) string
end
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded.");

diameterOrder = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"];
phi2Deg = (0:10:60)';
kernelHalfwidthDeg = 5;
gridResolutionDeg = 1;
displayPhi1Deg = 0:gridResolutionDeg:90;
displayPhiDeg = 0:gridResolutionDeg:90;
completePhi1Deg = 0:gridResolutionDeg:359;
completePhiDeg = displayPhiDeg;
[displayPhi1Grid,displayPhiGrid] = meshgrid( ...
  displayPhi1Deg,displayPhiDeg);
[completePhi1Grid,completePhiGrid] = meshgrid( ...
  completePhi1Deg,completePhiDeg);

catalog = comprehensive_ebsd_catalog(scanRoot);
catalog = catalog(catalog.variant == "raw",:);
[isRegistered,catalogOrder] = ismember(diameterOrder,catalog.sample);
assert(all(isRegistered) && numel(unique(catalogOrder)) == 6);
catalog = catalog(catalogOrder,:);
assert(height(catalog) == 6 && isequal(catalog.sample,diameterOrder));

kernel = SO3DeLaValleePoussinKernel( ...
  "halfwidth",kernelHalfwidthDeg * degree);
sampleCount = height(catalog);
sectionCount = numel(phi2Deg);
recordCount = sampleCount * sectionCount;
odfs = cell(sampleCount,1);
squareSectionValues = cell(sampleCount,sectionCount);
validCounts = zeros(sampleCount,1);
displayPeakMrd = zeros(recordCount,1);
displayPhi1PeakDeg = zeros(recordCount,1);
displayPhiPeakDeg = zeros(recordCount,1);
completePeakMrd = zeros(recordCount,1);
completePhi1PeakDeg = zeros(recordCount,1);
completePhiPeakDeg = zeros(recordCount,1);
cAxisAnglesDeg = zeros(recordCount,3);
orientationNote = strings(recordCount,1);
axisNames = ["AD","TD/RD","ND"];

for sampleIndex = 1:sampleCount
  [ebsdFull,~] = load_comprehensive_ebsd_scan(catalog(sampleIndex,:));
  orientations = ebsdFull("Ti-Hex").orientations;
  assert(~isempty(orientations));
  assert(string(orientations.CS.pointGroup) == "6/mmm");
  assert(string(orientations.CS.properGroup.pointGroup) == "622");
  assert(string(orientations.SS.pointGroup) == "1");
  [maxPhi1,maxPhi,maxPhi2] = fundamentalRegionEuler( ...
    orientations.CS,orientations.SS);
  assert(isequal([maxPhi1,maxPhi,maxPhi2],[360 90 60] * degree));
  validCounts(sampleIndex) = numel(orientations);

  rbfOdf = calcDensity(orientations,"kernel",kernel, ...
    "weights",ones(numel(orientations),1),"silent");
  rbfOdf = normalize_positive_mean_density(rbfOdf);
  odf = SO3FunHarmonic(rbfOdf,"bandwidth",kernel.bandwidth);
  odf = normalize_positive_mean_density(odf);
  assert(string(odf.CS.pointGroup) == "6/mmm");
  assert(string(odf.CS.properGroup.pointGroup) == "622");
  assert(string(odf.SS.pointGroup) == "1");
  odfs{sampleIndex} = odf;

  for sectionIndex = 1:sectionCount
    rowIndex = (sampleIndex - 1) * sectionCount + sectionIndex;
    sectionPhi2 = phi2Deg(sectionIndex);
    displayOrientations = orientation.byEuler( ...
      displayPhi1Grid(:) * degree,displayPhiGrid(:) * degree, ...
      sectionPhi2 * degree,"Bunge",odf.CS,odf.SS);
    displayValues = real(eval(odf,displayOrientations));
    assert(all(isfinite(displayValues)));
    displayValues = reshape(displayValues,size(displayPhi1Grid));
    squareSectionValues{sampleIndex,sectionIndex} = displayValues;
    [displayPeakMrd(rowIndex),displayIndex] = max(displayValues(:));
    displayPhi1PeakDeg(rowIndex) = displayPhi1Grid(displayIndex);
    displayPhiPeakDeg(rowIndex) = displayPhiGrid(displayIndex);

    completeOrientations = orientation.byEuler( ...
      completePhi1Grid(:) * degree,completePhiGrid(:) * degree, ...
      sectionPhi2 * degree,"Bunge",odf.CS,odf.SS);
    completeValues = real(eval(odf,completeOrientations));
    assert(all(isfinite(completeValues)));
    [completePeakMrd(rowIndex),completeIndex] = max(completeValues);
    completePhi1PeakDeg(rowIndex) = completePhi1Grid(completeIndex);
    completePhiPeakDeg(rowIndex) = completePhiGrid(completeIndex);

    peakOrientation = orientation.byEuler( ...
      displayPhi1PeakDeg(rowIndex) * degree, ...
      displayPhiPeakDeg(rowIndex) * degree,sectionPhi2 * degree, ...
      "Bunge",odf.CS,odf.SS);
    cDirection = peakOrientation * Miller(0,0,0,1,odf.CS);
    cAxisAnglesDeg(rowIndex,:) = [ ...
      angle(cDirection,vector3d.X,"antipodal"), ...
      angle(cDirection,vector3d.Y,"antipodal"), ...
      angle(cDirection,vector3d.Z,"antipodal")] / degree;
    [minimumAngle,minimumIndex] = min(cAxisAnglesDeg(rowIndex,:));
    orientationNote(rowIndex) = sprintf( ...
      "c-axis closest to %s (%.1f deg)", ...
      axisNames(minimumIndex),minimumAngle);
  end
  clear ebsdFull orientations rbfOdf odf
end

assert(all(isfinite(displayPeakMrd) & displayPeakMrd > 0));
assert(all(isfinite(completePeakMrd) & completePeakMrd > 0));
assert(all(completePeakMrd + 1e-10 >= displayPeakMrd));
globalMaximumMrd = max(displayPeakMrd);

peak_id = compose("P%02d",(1:recordCount)');
sample = repelem(catalog.sample,sectionCount);
diameter_mm = repelem(catalog.diameter_mm,sectionCount);
cold_reduction_percent = repelem( ...
  catalog.cold_reduction_percent,sectionCount);
input_path = repelem(catalog.input_path,sectionCount);
valid_ti_hex_orientation_count = repelem(validCounts,sectionCount);
phi2_deg = repmat(phi2Deg,sampleCount,1);
display_peak_mrd = displayPeakMrd;
display_phi1_deg = displayPhi1PeakDeg;
display_Phi_deg = displayPhiPeakDeg;
display_phi2_deg = phi2_deg;
c_axis_to_ad_deg = cAxisAnglesDeg(:,1);
c_axis_to_td_rd_deg = cAxisAnglesDeg(:,2);
c_axis_to_nd_deg = cAxisAnglesDeg(:,3);
orientation_note = orientationNote;
complete_peak_mrd = completePeakMrd;
complete_phi1_deg = completePhi1PeakDeg;
complete_Phi_deg = completePhiPeakDeg;
complete_phi2_deg = phi2_deg;
complete_peak_inside_display = complete_phi1_deg <= 90;
complete_minus_display_mrd = complete_peak_mrd - display_peak_mrd;
global_color_limit_max_mrd = repmat(globalMaximumMrd,recordCount,1);
peakSummary = table(peak_id,sample,diameter_mm, ...
  cold_reduction_percent,input_path,valid_ti_hex_orientation_count, ...
  phi2_deg,display_peak_mrd,display_phi1_deg,display_Phi_deg, ...
  display_phi2_deg,c_axis_to_ad_deg,c_axis_to_td_rd_deg, ...
  c_axis_to_nd_deg,orientation_note,complete_peak_mrd, ...
  complete_phi1_deg,complete_Phi_deg,complete_phi2_deg, ...
  complete_peak_inside_display,complete_minus_display_mrd, ...
  global_color_limit_max_mrd);

data_variant = "raw";
phase = "Ti-Hex";
weighting = "pixel";
crystal_symmetry = "6/mmm";
rotational_group = "622";
specimen_symmetry = "1";
x_axis = "AD";
y_axis = "TD/RD";
z_axis = "ND";
calculation_method = "harmonic series expansion (SO3FunHarmonic)";
kernel_type = "De la Vallee Poussin";
kernel_halfwidth_deg = kernelHalfwidthDeg;
grid_resolution_deg = gridResolutionDeg;
phi2_sections_deg = "0,10,20,30,40,50,60";
display_phi1_range_deg = "0:1:90";
display_Phi_range_deg = "0:1:90";
complete_phi1_range_deg = "0:1:359";
complete_Phi_range_deg = "0:1:90";
normalization = "unit mean MRD";
tie_breaking_rule = "lowest Phi, then lowest phi1";
color_limit_definition = "maximum of 42 display-domain peaks";
parameterSummary = table(data_variant,phase,weighting, ...
  crystal_symmetry,rotational_group,specimen_symmetry,x_axis,y_axis, ...
  z_axis,calculation_method,kernel_type,kernel_halfwidth_deg, ...
  grid_resolution_deg,phi2_sections_deg,display_phi1_range_deg, ...
  display_Phi_range_deg,complete_phi1_range_deg, ...
  complete_Phi_range_deg,normalization,tie_breaking_rule, ...
  color_limit_definition);

assert(height(peakSummary) == 42 && height(parameterSummary) == 1);
assert(isequal(peakSummary.peak_id,compose("P%02d",(1:42)')));
end
```

- [ ] **Step 4: Run the calculation test and verify it passes**

```powershell
$batch = "addpath(fullfile(pwd,'tools','mtex')); startup_project_mtex; scanRoot=string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')); test_calculate_square_odf_peak_data(scanRoot);"
Invoke-ProjectMatlabBatch -Batch $batch -LogStem 'square_odf_calculation_green'
```

Expected: exit code 0 and
`test_calculate_square_odf_peak_data passed`.

- [ ] **Step 5: Commit the tested calculation layer**

```powershell
git add -- 'tools/mtex/calculate_square_odf_peak_data.m' 'tools/mtex/test_calculate_square_odf_peak_data.m'
git commit -m "feat: calculate square ODF peak orientations"
```

---

### Task 2: Render and export the square peak matrix

**Files:**
- Create: `tools/mtex/test_generate_odf_square_peak_matrix.m`
- Create: `tools/mtex/generate_odf_square_peak_matrix.m`

**Interfaces:**
- Consumes: all five outputs of `calculate_square_odf_peak_data(scanRoot)`.
- Produces: `[peakSummary,parameterSummary] = generate_odf_square_peak_matrix(scanRoot,outputRoot)`.
- Writes exactly `odf_square_peak_matrix.png`, `odf_square_peak_matrix.pdf`, `odf_square_peak_summary.csv`, and `odf_square_peak_parameters.csv`.

- [ ] **Step 1: Write the failing generator test**

Create `tools/mtex/test_generate_odf_square_peak_matrix.m`:

```matlab
function test_generate_odf_square_peak_matrix(scanRoot)
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded for this test.");
outputRoot = string(tempname);
mkdir(outputRoot);
cleanupOutput = onCleanup(@() remove_test_output(outputRoot));
[peaks,parameters] = ...
  generate_odf_square_peak_matrix(scanRoot,outputRoot);

expectedPeakColumns = ["peak_id","sample","diameter_mm", ...
  "cold_reduction_percent","input_path", ...
  "valid_ti_hex_orientation_count","phi2_deg","display_peak_mrd", ...
  "display_phi1_deg","display_Phi_deg","display_phi2_deg", ...
  "c_axis_to_ad_deg","c_axis_to_td_rd_deg","c_axis_to_nd_deg", ...
  "orientation_note","complete_peak_mrd","complete_phi1_deg", ...
  "complete_Phi_deg","complete_phi2_deg", ...
  "complete_peak_inside_display","complete_minus_display_mrd", ...
  "global_color_limit_max_mrd"];
expectedParameterColumns = ["data_variant","phase","weighting", ...
  "crystal_symmetry","rotational_group","specimen_symmetry", ...
  "x_axis","y_axis","z_axis","calculation_method","kernel_type", ...
  "kernel_halfwidth_deg","grid_resolution_deg","phi2_sections_deg", ...
  "display_phi1_range_deg","display_Phi_range_deg", ...
  "complete_phi1_range_deg","complete_Phi_range_deg", ...
  "normalization","tie_breaking_rule","color_limit_definition"];
assert(isequal(string(peaks.Properties.VariableNames), ...
  expectedPeakColumns));
assert(isequal(string(parameters.Properties.VariableNames), ...
  expectedParameterColumns));
assert(height(peaks) == 42 && height(parameters) == 1);

peakCsv = fullfile(outputRoot,"odf_square_peak_summary.csv");
parameterCsv = fullfile(outputRoot,"odf_square_peak_parameters.csv");
pngPath = fullfile(outputRoot,"odf_square_peak_matrix.png");
pdfPath = fullfile(outputRoot,"odf_square_peak_matrix.pdf");
assert(isfile(peakCsv) && isfile(parameterCsv));
assert(isfile(pngPath) && isfile(pdfPath));
assert(height(readtable(peakCsv,"TextType","string")) == 42);
assert(height(readtable(parameterCsv,"TextType","string")) == 1);

imageInfo = imfinfo(pngPath);
assert(imageInfo.Width > imageInfo.Height);
assert(imageInfo.Width >= 6000 && imageInfo.Height >= 4500);
[xDpi,yDpi] = image_resolution_dpi(imageInfo);
assert(abs(xDpi - 600) < 1 && abs(yDpi - 600) < 1);
imageData = imread(pngPath);
assert(mean(any(imageData < 250,3),"all") >= 0.15);

pdfInfo = dir(pdfPath);
assert(pdfInfo.bytes > 0);
pdfText = fileread(pdfPath);
assert(startsWith(pdfText,"%PDF-"));
imageTokens = regexp(pdfText, ...
  '/Subtype /Image\s+/Width ([0-9]+)\s+/Height ([0-9]+)', ...
  "tokens");
mediaBox = regexp(pdfText, ...
  '/MediaBox \[0 0 ([0-9.]+) ([0-9.]+)\]',"tokens","once");
assert(~isempty(imageTokens) && ~isempty(mediaBox));
pixelArea = sum(cellfun(@(token) ...
  str2double(token{1}) * str2double(token{2}),imageTokens));
pageArea = (str2double(mediaBox{1}) / 72) * ...
  (str2double(mediaBox{2}) / 72);
assert(sqrt(pixelArea / pageArea) >= 590);
clear cleanupOutput
fprintf("test_generate_odf_square_peak_matrix passed\n");
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

- [ ] **Step 2: Run the generator test and observe the expected failure**

```powershell
$batch = "addpath(fullfile(pwd,'tools','mtex')); startup_project_mtex; scanRoot=string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')); test_generate_odf_square_peak_matrix(scanRoot);"
Invoke-ProjectMatlabBatch -Batch $batch -LogStem 'square_odf_generator_red'
```

Expected: MATLAB exits nonzero because
`generate_odf_square_peak_matrix` is undefined.

- [ ] **Step 3: Implement the generator and renderer**

Create `tools/mtex/generate_odf_square_peak_matrix.m`:

```matlab
function [peakSummary,parameterSummary] = ...
  generate_odf_square_peak_matrix(scanRoot,outputRoot)
%GENERATE_ODF_SQUARE_PEAK_MATRIX Render square alpha-Ti ODF sections.

arguments
  scanRoot (1,1) string
  outputRoot (1,1) string
end
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded.");
if ~isfolder(outputRoot)
  mkdir(outputRoot);
end
[peakSummary,parameterSummary,~,catalog,squareValues] = ...
  calculate_square_odf_peak_data(scanRoot);
writetable(peakSummary, ...
  fullfile(outputRoot,"odf_square_peak_summary.csv"));
writetable(parameterSummary, ...
  fullfile(outputRoot,"odf_square_peak_parameters.csv"));
pngPath = fullfile(outputRoot,"odf_square_peak_matrix.png");
pdfPath = fullfile(outputRoot,"odf_square_peak_matrix.pdf");
render_square_peak_matrix(squareValues,peakSummary,catalog,pngPath,pdfPath);
assert(isfile(pngPath) && isfile(pdfPath));
end

function render_square_peak_matrix(squareValues,peaks,catalog, ...
    pngPath,pdfPath)
sampleCount = 6;
sectionCount = 7;
assert(isequal(size(squareValues),[sampleCount sectionCount]));
assert(height(peaks) == sampleCount * sectionCount);
globalLimit = peaks.global_color_limit_max_mrd(1);
assert(all(peaks.global_color_limit_max_mrd == globalLimit));
contourLevels = linspace(0,globalLimit,13);
phi1Deg = 0:90;
PhiDeg = 0:90;

figureHandle = figure("Visible","off","Color","white", ...
  "Units","inches","Position",[0.25 0.25 15.2 12.5]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle,sampleCount,sectionCount, ...
  "TileSpacing","compact","Padding","compact");
axesHandles = gobjects(sampleCount,sectionCount);
for sampleIndex = 1:sampleCount
  for sectionIndex = 1:sectionCount
    rowIndex = (sampleIndex - 1) * sectionCount + sectionIndex;
    axesHandle = nexttile(layout,rowIndex);
    axesHandles(sampleIndex,sectionIndex) = axesHandle;
    values = max(real(squareValues{sampleIndex,sectionIndex}),0);
    assert(isequal(size(values),[91 91]) && all(isfinite(values),"all"));
    contourf(axesHandle,phi1Deg,PhiDeg,values,contourLevels, ...
      "LineColor",[0.20 0.20 0.20],"LineWidth",0.25);
    hold(axesHandle,"on");
    plot(axesHandle,peaks.display_phi1_deg(rowIndex), ...
      peaks.display_Phi_deg(rowIndex),"wo", ...
      "MarkerSize",5,"LineWidth",1.2);
    plot(axesHandle,peaks.display_phi1_deg(rowIndex), ...
      peaks.display_Phi_deg(rowIndex),"k+", ...
      "MarkerSize",5,"LineWidth",0.8);
    text(axesHandle,0.03,0.97,sprintf("%s | %.2f", ...
      peaks.peak_id(rowIndex),peaks.display_peak_mrd(rowIndex)), ...
      "Units","normalized","HorizontalAlignment","left", ...
      "VerticalAlignment","top","FontSize",5.5, ...
      "FontWeight","bold","Color","black", ...
      "BackgroundColor","white","Margin",0.5);
    axis(axesHandle,"xy");
    axis(axesHandle,"square");
    xlim(axesHandle,[0 90]);
    ylim(axesHandle,[0 90]);
    clim(axesHandle,[0 globalLimit]);
    xticks(axesHandle,[0 45 90]);
    yticks(axesHandle,[0 45 90]);
    axesHandle.FontSize = 5.5;
    axesHandle.LineWidth = 0.5;
    if sampleIndex < sampleCount
      axesHandle.XTickLabel = [];
    else
      xlabel(axesHandle,"\phi_1 (deg)","FontSize",6.5);
    end
    if sectionIndex > 1
      axesHandle.YTickLabel = [];
    else
      ylabel(axesHandle,sprintf( ...
        "%.2f mm | %.2f%%\\n\\Phi (deg)", ...
        catalog.diameter_mm(sampleIndex), ...
        catalog.cold_reduction_percent(sampleIndex)), ...
        "FontSize",6.5);
    end
    if sampleIndex == 1
      title(axesHandle,sprintf("\\phi_2 = %d^\\circ", ...
        peaks.phi2_deg(rowIndex)),"FontSize",7.5);
    end
  end
end
colormap(figureHandle,parula(256));
colorbarHandle = colorbar(axesHandles(end,end));
colorbarHandle.Layout.Tile = "east";
colorbarHandle.Label.String = "ODF intensity (MRD)";
colorbarHandle.FontSize = 7;
title(layout, ...
  "Square ODF sections (SS = 1; displayed \phi_1 = 0-90^\circ)", ...
  "FontSize",10,"FontWeight","bold");
xlabel(layout, ...
  "Markers and labels P01-P42 denote maxima within the displayed square region", ...
  "FontSize",7);
drawnow;
exportgraphics(figureHandle,char(pngPath), ...
  "Resolution",600,"BackgroundColor","white");
exportgraphics(figureHandle,char(pdfPath), ...
  "ContentType","image","Resolution",600,"BackgroundColor","white");
clear cleanupFigure
end
```

- [ ] **Step 4: Run the generator test and verify it passes**

```powershell
$batch = "addpath(fullfile(pwd,'tools','mtex')); startup_project_mtex; scanRoot=string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')); test_generate_odf_square_peak_matrix(scanRoot);"
Invoke-ProjectMatlabBatch -Batch $batch -LogStem 'square_odf_generator_green'
```

Expected: exit code 0 and
`test_generate_odf_square_peak_matrix passed`.

- [ ] **Step 5: Run static checks and the focused regression suite**

```powershell
$batch = "addpath(fullfile(pwd,'tools','mtex')); startup_project_mtex; scanRoot=string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')); test_comprehensive_ebsd_contract; test_calculate_square_odf_peak_data(scanRoot); test_generate_odf_square_peak_matrix(scanRoot); files={'tools/mtex/calculate_square_odf_peak_data.m','tools/mtex/test_calculate_square_odf_peak_data.m','tools/mtex/generate_odf_square_peak_matrix.m','tools/mtex/test_generate_odf_square_peak_matrix.m'}; for k=1:numel(files), issues=checkcode(files{k},'-id'); assert(isempty(issues),files{k}); end; disp('SQUARE_ODF_FOCUSED_SUITE_OK');"
Invoke-ProjectMatlabBatch -Batch $batch -LogStem 'square_odf_focused_suite'
```

Expected: exit code 0 and `SQUARE_ODF_FOCUSED_SUITE_OK`.

- [ ] **Step 6: Commit the tested generator**

```powershell
git add -- 'tools/mtex/generate_odf_square_peak_matrix.m' 'tools/mtex/test_generate_odf_square_peak_matrix.m'
git commit -m "feat: render square ODF peak matrix"
```

---

### Task 3: Generate, inspect, and record the production outputs

**Files:**
- Create: `results/mtex_odf_square_peak_matrix/odf_square_peak_matrix.png`
- Create: `results/mtex_odf_square_peak_matrix/odf_square_peak_matrix.pdf`
- Create: `results/mtex_odf_square_peak_matrix/odf_square_peak_summary.csv`
- Create: `results/mtex_odf_square_peak_matrix/odf_square_peak_parameters.csv`

**Interfaces:**
- Consumes: the tested `generate_odf_square_peak_matrix` entry point.
- Produces: the four reviewable publication/audit artifacts in the approved result directory.

- [ ] **Step 1: Generate the production outputs in a hidden MATLAB batch**

```powershell
$batch = "addpath(fullfile(pwd,'tools','mtex')); startup_project_mtex; scanRoot=string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')); outputRoot=string(fullfile(pwd,'results','mtex_odf_square_peak_matrix')); generate_odf_square_peak_matrix(scanRoot,outputRoot); disp('SQUARE_ODF_PRODUCTION_OK');"
Invoke-ProjectMatlabBatch -Batch $batch -LogStem 'square_odf_production'
```

Expected: exit code 0 and `SQUARE_ODF_PRODUCTION_OK`.

- [ ] **Step 2: Audit the CSVs and file inventory**

```powershell
$resultRoot = 'results\mtex_odf_square_peak_matrix'
$peaks = Import-Csv -LiteralPath (Join-Path $resultRoot 'odf_square_peak_summary.csv')
$parameters = Import-Csv -LiteralPath (Join-Path $resultRoot 'odf_square_peak_parameters.csv')
if ($peaks.Count -ne 42) { throw "Expected 42 peak rows." }
if ($parameters.Count -ne 1) { throw "Expected one parameter row." }
if (($peaks.peak_id -join ',') -ne
    ((1..42 | ForEach-Object { 'P{0:D2}' -f $_ }) -join ',')) {
  throw "Peak identifiers are not P01-P42."
}
if (($peaks.phi2_deg | Sort-Object -Unique) -join ',' -ne
    '0,10,20,30,40,50,60') {
  throw "Unexpected phi2 sections."
}
if (($peaks | Where-Object {
      [double]$_.complete_peak_mrd + 1e-10 -lt
      [double]$_.display_peak_mrd
    }).Count -ne 0) {
  throw "A display-domain peak exceeds its complete-section peak."
}
Get-ChildItem -LiteralPath $resultRoot -File |
  Select-Object Name,Length,LastWriteTime
```

Expected: four non-empty files and no thrown assertion.

- [ ] **Step 3: Visually inspect the PNG**

Open
`results/mtex_odf_square_peak_matrix/odf_square_peak_matrix.png`
and verify:

- six diameter rows and seven `phi2 = 0:10:60 deg` columns;
- every contour panel is square and spans `0--90 deg` on both plotted axes;
- `P01` through `P42` are present in row-major order;
- every panel contains a visible peak marker and MRD value;
- row labels contain diameter and cold reduction;
- the shared colorbar is readable and no plotted peak is clipped;
- there are no blank panels, contour holes, overlapping labels, or cropped
  outer margins.

If a visual defect is found, add a focused assertion where practical, change
only the renderer, rerun Task 2 Step 5, regenerate, and repeat this inspection.

- [ ] **Step 4: Run the final verification batch against production outputs**

```powershell
$batch = "addpath(fullfile(pwd,'tools','mtex')); startup_project_mtex; scanRoot=string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')); test_comprehensive_ebsd_contract; test_calculate_square_odf_peak_data(scanRoot); test_generate_odf_square_peak_matrix(scanRoot); outputRoot=string(fullfile(pwd,'results','mtex_odf_square_peak_matrix')); assert(isfile(fullfile(outputRoot,'odf_square_peak_matrix.png'))); assert(isfile(fullfile(outputRoot,'odf_square_peak_matrix.pdf'))); peaks=readtable(fullfile(outputRoot,'odf_square_peak_summary.csv'),'TextType','string'); parameters=readtable(fullfile(outputRoot,'odf_square_peak_parameters.csv'),'TextType','string'); assert(height(peaks)==42 && height(parameters)==1); assert(all(peaks.complete_peak_mrd+1e-10>=peaks.display_peak_mrd)); disp('SQUARE_ODF_FINAL_VERIFICATION_OK');"
Invoke-ProjectMatlabBatch -Batch $batch -LogStem 'square_odf_final_verification'
```

Expected: exit code 0 and `SQUARE_ODF_FINAL_VERIFICATION_OK`.

- [ ] **Step 5: Commit only the four production artifacts**

```powershell
git add -- `
  'results/mtex_odf_square_peak_matrix/odf_square_peak_matrix.png' `
  'results/mtex_odf_square_peak_matrix/odf_square_peak_matrix.pdf' `
  'results/mtex_odf_square_peak_matrix/odf_square_peak_summary.csv' `
  'results/mtex_odf_square_peak_matrix/odf_square_peak_parameters.csv'
git diff --cached --check
git commit -m "results: add square ODF peak matrix"
git status --short
```

Expected: the result commit contains exactly four files; unrelated worktree
changes, if any, remain untouched.
