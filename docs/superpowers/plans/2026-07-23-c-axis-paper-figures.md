# C-Axis Paper Figures Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate separate raw and denoised 2×3 `{0001}` c-axis pole-figure montages and separate raw and denoised c-axis–AD mean-orientation trend figures.

**Architecture:** Add one pure table-preparation function and one focused MTEX figure generator. The preparation function validates and sorts the existing pixel-weighted texture summary; the generator rebuilds only the 12 pixel-weighted ODFs needed for `{0001}` pole figures, computes one shared MRD limit, and writes four new figures without touching existing artifacts or source CTF files.

**Tech Stack:** MATLAB R2025a, MTEX 6.1.1, existing comprehensive EBSD catalog/load helpers, MATLAB `tiledlayout`, `plotPDF`, and `exportgraphics`.

## Global Constraints

- Use the six states in this exact order: 0%, 14.31%, 26.04%, 36%, 43.75%, 48.98%.
- Produce separate raw and denoised figures; never treat denoised data as an independent repeat.
- Use pixel-weighted ODFs with the registered 5° de la Vallée Poussin kernel and 5° evaluation grid.
- Use `{0001}`, upper-hemisphere equal-area antipodal pole figures with AD horizontal.
- Use one shared MRD range across all 12 pole figures.
- Use identical axis limits for the separate raw and denoised mean-orientation plots.
- Show the mean c-axis–AD acute angle and its P10–P90 range.
- Do not overwrite existing figures or any CTF input.
- Write all four outputs to `results/mtex_ebsd_comprehensive/05_texture/` at 300 dpi.

---

### Task 1: Validate and prepare mean-orientation rows

**Files:**
- Create: `tools/mtex/prepare_c_axis_figure_rows.m`
- Create: `tools/mtex/test_c_axis_paper_figures.m`

**Interfaces:**
- Consumes: `textureSummary table` with the registered texture-summary columns.
- Produces: `prepared struct` with fields `raw`, `denoised`, `cold_reduction_percent`, and `common_y_limits_deg`.

- [ ] **Step 1: Write the failing pure-table test**

Add this test entry point and fixture to `test_c_axis_paper_figures.m`:

```matlab
function test_c_axis_paper_figures(scanRoot, outputRoot)
arguments
  scanRoot (1,1) string = ""
  outputRoot (1,1) string = ""
end

test_prepare_rows();
if scanRoot ~= ""
  assert(outputRoot ~= "");
  test_formal_generation(scanRoot, outputRoot);
end
fprintf("test_c_axis_paper_figures passed\n");
end

function test_prepare_rows()
sample = repelem(["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"], 4);
diameter_mm = repelem([7;6.48;6.02;5.6;5.25;5], 4);
cold_reduction_percent = repelem([0;14.31;26.04;36;43.75;48.98], 4);
variant = repmat(["raw";"raw";"denoised";"denoised"], 6, 1);
weighting = repmat(["pixel_weighted";"area_weighted_grain_mean"; ...
  "pixel_weighted";"area_weighted_grain_mean"], 6, 1);
c_axis_ad_mean_deg = repelem((79:84)',4) + ...
  repmat([0;20;0.1;20.1],6,1);
c_axis_ad_p10_deg = c_axis_ad_mean_deg - 10;
c_axis_ad_p90_deg = c_axis_ad_mean_deg + 8;
summary = table(sample,diameter_mm,cold_reduction_percent,variant, ...
  weighting,c_axis_ad_mean_deg,c_axis_ad_p10_deg,c_axis_ad_p90_deg);
prepared = prepare_c_axis_figure_rows(summary);
assert(height(prepared.raw) == 6 && height(prepared.denoised) == 6);
assert(all(prepared.raw.weighting == "pixel_weighted"));
assert(isequal(prepared.cold_reduction_percent, ...
  [0;14.31;26.04;36;43.75;48.98]));
assert(all(abs(prepared.denoised.c_axis_ad_mean_deg - ...
  prepared.raw.c_axis_ad_mean_deg - 0.1) < 1e-12));
assert(prepared.common_y_limits_deg(1) <= ...
  min([prepared.raw.c_axis_ad_p10_deg; ...
  prepared.denoised.c_axis_ad_p10_deg]));
assert(prepared.common_y_limits_deg(2) >= ...
  max([prepared.raw.c_axis_ad_p90_deg; ...
  prepared.denoised.c_axis_ad_p90_deg]));
end
```

- [ ] **Step 2: Run the test to verify RED**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath('C:/Users/22069/Documents/MATLAB/mtex-6.1.1'); startup_mtex; addpath('tools/mtex'); test_c_axis_paper_figures"
```

Expected: FAIL because `prepare_c_axis_figure_rows` is undefined.

- [ ] **Step 3: Implement the pure preparation function**

Create `prepare_c_axis_figure_rows.m`:

```matlab
function prepared = prepare_c_axis_figure_rows(textureSummary)
required = ["sample","diameter_mm","cold_reduction_percent", ...
  "variant","weighting","c_axis_ad_mean_deg", ...
  "c_axis_ad_p10_deg","c_axis_ad_p90_deg"];
assert(istable(textureSummary));
assert(all(ismember(required, ...
  string(textureSummary.Properties.VariableNames))));
textureSummary.sample = string(textureSummary.sample);
textureSummary.variant = string(textureSummary.variant);
textureSummary.weighting = string(textureSummary.weighting);
rows = textureSummary.weighting == "pixel_weighted";
selected = textureSummary(rows, cellstr(required));
assert(height(selected) == 12);
key = selected.sample + "|" + selected.variant;
assert(numel(unique(key)) == 12);
expectedReduction = [0;14.31;26.04;36;43.75;48.98];
prepared = struct();
for variantName = ["raw","denoised"]
  variantRows = selected(selected.variant == variantName,:);
  variantRows = sortrows(variantRows,"cold_reduction_percent");
  assert(height(variantRows) == 6);
  assert(max(abs(variantRows.cold_reduction_percent - ...
    expectedReduction)) < 1e-10);
  assert(all(isfinite(variantRows{:,6:8}),"all"));
  prepared.(variantName) = variantRows;
end
prepared.cold_reduction_percent = expectedReduction;
lower = min([prepared.raw.c_axis_ad_p10_deg; ...
  prepared.denoised.c_axis_ad_p10_deg]);
upper = max([prepared.raw.c_axis_ad_p90_deg; ...
  prepared.denoised.c_axis_ad_p90_deg]);
padding = max(1,0.05 * (upper - lower));
prepared.common_y_limits_deg = ...
  [floor(lower-padding),ceil(upper+padding)];
end
```

- [ ] **Step 4: Run the pure test to verify GREEN**

Run the command from Step 2.

Expected: `test_c_axis_paper_figures passed`.

- [ ] **Step 5: Commit Task 1**

```powershell
git add -- tools/mtex/prepare_c_axis_figure_rows.m tools/mtex/test_c_axis_paper_figures.m
git commit -m "test: define c-axis paper figure data contract"
```

---

### Task 2: Generate the four publication figures

**Files:**
- Create: `tools/mtex/generate_c_axis_paper_figures.m`
- Modify: `tools/mtex/test_c_axis_paper_figures.m`

**Interfaces:**
- Consumes: `scanRoot string`, `outputRoot string`, the existing `05_texture/texture_summary.csv`, and the registered catalog.
- Produces: `metadata table` with one row per output and the common MRD/y-axis limits; writes the four approved PNG files.

- [ ] **Step 1: Extend the test with a formal-output assertion**

Add:

```matlab
function test_formal_generation(scanRoot, outputRoot)
catalog = comprehensive_ebsd_catalog(scanRoot);
[beforeBytes,beforeTimes] = input_stats(catalog.input_path);
metadata = generate_c_axis_paper_figures(scanRoot,outputRoot);
expected = ["c_axis_pole_figures_raw.png"; ...
  "c_axis_pole_figures_denoised.png"; ...
  "c_axis_mean_orientation_raw.png"; ...
  "c_axis_mean_orientation_denoised.png"];
assert(isequal(metadata.output_name,expected));
assert(height(metadata) == 4);
assert(numel(unique(metadata.pole_color_max_mrd(1:2))) == 1);
assert(numel(unique(metadata.y_min_deg(3:4))) == 1);
assert(numel(unique(metadata.y_max_deg(3:4))) == 1);
textureDir = fullfile(outputRoot,"05_texture");
for fileName = expected'
  info = dir(fullfile(textureDir,fileName));
  assert(isscalar(info) && info.bytes > 0);
end
[afterBytes,afterTimes] = input_stats(catalog.input_path);
assert(isequal(beforeBytes,afterBytes));
assert(isequal(beforeTimes,afterTimes));
end

function [bytes,times] = input_stats(paths)
bytes = zeros(numel(paths),1);
times = zeros(numel(paths),1);
for pathIndex = 1:numel(paths)
  info = dir(paths(pathIndex));
  assert(isscalar(info));
  bytes(pathIndex) = info.bytes;
  times(pathIndex) = info.datenum;
end
end
```

- [ ] **Step 2: Run the full test to verify RED**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath('C:/Users/22069/Documents/MATLAB/mtex-6.1.1'); startup_mtex; addpath('tools/mtex'); test_c_axis_paper_figures('data/ebsd_kpl_250221_7_df/scans','results/mtex_ebsd_comprehensive')"
```

Expected: FAIL because `generate_c_axis_paper_figures` is undefined.

- [ ] **Step 3: Implement ODF rebuilding and shared-scale calculation**

Create `generate_c_axis_paper_figures.m` with:

```matlab
function metadata = generate_c_axis_paper_figures(scanRoot,outputRoot)
arguments
  scanRoot (1,1) string
  outputRoot (1,1) string
end
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")));
textureDir = fullfile(outputRoot,"05_texture");
summaryPath = fullfile(textureDir,"texture_summary.csv");
assert(isfile(summaryPath));
summary = readtable(summaryPath,"TextType","string", ...
  "VariableNamingRule","preserve");
prepared = prepare_c_axis_figure_rows(summary);
catalog = comprehensive_ebsd_catalog(scanRoot);
assert(height(catalog) == 12);
parameters = comprehensive_ebsd_output_contract().parameters;
kernel = SO3DeLaValleePoussinKernel("halfwidth", ...
  parameters.texture_kernel_halfwidth_deg * degree);
odfs = cell(12,1);
sharedMaximum = 0;
for scanIndex = 1:12
  [ebsdFull,~] = load_comprehensive_ebsd_scan(catalog(scanIndex,:));
  tiEbsd = ebsdFull("Ti-Hex");
  rbfOdf = calcDensity(tiEbsd.orientations,"kernel",kernel, ...
    "weights",ones(length(tiEbsd),1),"silent");
  rbfOdf = rbfOdf / double(mean(rbfOdf));
  odfs{scanIndex} = SO3FunHarmonic(rbfOdf, ...
    "bandwidth",kernel.bandwidth);
  cAxis = Miller(0,0,0,1,tiEbsd.CS);
  density = calcPDF(odfs{scanIndex},cAxis,[],"antipodal");
  sharedMaximum = max(sharedMaximum,double(max(density, ...
    "resolution",parameters.texture_grid_resolution_deg*degree)));
  clear ebsdFull tiEbsd rbfOdf density
end
assert(isfinite(sharedMaximum) && sharedMaximum > 0);
```

Then call focused local rendering functions for the two variants and two trend tables, and return a four-row metadata table:

```matlab
output_name = ["c_axis_pole_figures_raw.png"; ...
  "c_axis_pole_figures_denoised.png"; ...
  "c_axis_mean_orientation_raw.png"; ...
  "c_axis_mean_orientation_denoised.png"];
variant = ["raw";"denoised";"raw";"denoised"];
figure_kind = ["pole_montage";"pole_montage"; ...
  "mean_orientation";"mean_orientation"];
pole_color_max_mrd = [sharedMaximum;sharedMaximum;NaN;NaN];
y_min_deg = [NaN;NaN;repmat(prepared.common_y_limits_deg(1),2,1)];
y_max_deg = [NaN;NaN;repmat(prepared.common_y_limits_deg(2),2,1)];
metadata = table(output_name,variant,figure_kind, ...
  pole_color_max_mrd,y_min_deg,y_max_deg);
```

- [ ] **Step 4: Implement the 2×3 pole montage renderer**

Add a local function that:

```matlab
function render_pole_montage(odfs,catalog,variantName, ...
  colorMaximum,outputPath,parameters)
rows = find(catalog.variant == variantName);
assert(numel(rows) == 6);
[~,order] = sort(catalog.cold_reduction_percent(rows));
rows = rows(order);
figureHandle = figure("Visible","off","Color","white", ...
  "Position",[50 50 1800 1180]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle,2,3,"Padding","compact", ...
  "TileSpacing","compact");
for panelIndex = 1:6
  axesHandle = nexttile(layout);
  odf = odfs{rows(panelIndex)};
  cAxis = Miller(0,0,0,1,odf.CS);
  plotPDF(odf,cAxis,"antipodal","earea","contourf","silent", ...
    "parent",axesHandle, ...
    "resolution",parameters.texture_grid_resolution_deg*degree, ...
    "colorRange",[0 colorMaximum]);
  title(axesHandle,sprintf("%.2f%%", ...
    catalog.cold_reduction_percent(rows(panelIndex))), ...
    "FontWeight","normal");
end
colorbarHandle = colorbar(axesHandle,"eastoutside");
colorbarHandle.Label.String = "MRD";
title(layout,sprintf( ...
  "Alpha-Ti {0001} c-axis pole figures | %s | AD horizontal", ...
  variantName),"Interpreter","none");
exportgraphics(figureHandle,outputPath,"Resolution",300, ...
  "BackgroundColor","white");
clear cleanupFigure
end
```

Use one explicit `setMTEXpref("pfAnnotations",...)` definition matching X=AD, Y=TD-RD, Z=ND before rendering.

- [ ] **Step 5: Implement the separate mean-orientation renderer**

Add:

```matlab
function render_mean_trend(rows,variantName,yLimits,outputPath)
figureHandle = figure("Visible","off","Color","white", ...
  "Position",[100 100 1200 760]);
cleanupFigure = onCleanup(@() close(figureHandle));
axesHandle = axes(figureHandle);
hold(axesHandle,"on");
x = rows.cold_reduction_percent;
fill(axesHandle,[x;flipud(x)], ...
  [rows.c_axis_ad_p10_deg;flipud(rows.c_axis_ad_p90_deg)], ...
  [0.75 0.84 0.94],"EdgeColor","none","FaceAlpha",0.45, ...
  "DisplayName","P10-P90");
plot(axesHandle,x,rows.c_axis_ad_mean_deg,"-o", ...
  "LineWidth",1.8,"MarkerSize",7, ...
  "DisplayName","Mean c-axis-AD angle");
for pointIndex = 1:height(rows)
  text(axesHandle,x(pointIndex),rows.c_axis_ad_mean_deg(pointIndex), ...
    sprintf("  %.2f",rows.c_axis_ad_mean_deg(pointIndex)), ...
    "VerticalAlignment","bottom");
end
xlim(axesHandle,[0 50]);
ylim(axesHandle,yLimits);
xlabel(axesHandle,"Cold reduction (%)");
ylabel(axesHandle,"c-axis-AD acute angle (deg)");
title(axesHandle,sprintf( ...
  "Mean alpha-Ti c-axis orientation | %s | AD horizontal", ...
  variantName),"Interpreter","none","FontWeight","normal");
grid(axesHandle,"on");
legend(axesHandle,"Location","best");
exportgraphics(figureHandle,outputPath,"Resolution",300, ...
  "BackgroundColor","white");
clear cleanupFigure
end
```

- [ ] **Step 6: Run the formal generator test to verify GREEN**

Run the command from Step 2.

Expected: `test_c_axis_paper_figures passed`, with four non-empty PNG files.

- [ ] **Step 7: Commit Task 2**

```powershell
git add -- tools/mtex/generate_c_axis_paper_figures.m tools/mtex/test_c_axis_paper_figures.m
git commit -m "feat: add separate c-axis paper figures"
```

---

### Task 3: Verify scientific values and visual layout

**Files:**
- Modify: `results/mtex_ebsd_comprehensive/README.md`
- Verify: `results/mtex_ebsd_comprehensive/05_texture/c_axis_pole_figures_raw.png`
- Verify: `results/mtex_ebsd_comprehensive/05_texture/c_axis_pole_figures_denoised.png`
- Verify: `results/mtex_ebsd_comprehensive/05_texture/c_axis_mean_orientation_raw.png`
- Verify: `results/mtex_ebsd_comprehensive/05_texture/c_axis_mean_orientation_denoised.png`

**Interfaces:**
- Consumes: the four outputs and `texture_summary.csv`.
- Produces: visually verified manuscript-ready figures and README registration.

- [ ] **Step 1: Compare plotted means with the registered summary**

Run:

```powershell
$rows = Import-Csv 'results/mtex_ebsd_comprehensive/05_texture/texture_summary.csv' |
  Where-Object { $_.weighting -eq 'pixel_weighted' } |
  Sort-Object variant,{[double]$_.cold_reduction_percent}
$rows | Select-Object sample,variant,cold_reduction_percent,
  c_axis_ad_mean_deg,c_axis_ad_p10_deg,c_axis_ad_p90_deg |
  Format-Table -AutoSize
```

Expected: six finite ordered rows for each variant; the raw means remain approximately 78.6°–79.8°.

- [ ] **Step 2: Visually inspect all four PNG files**

Open each output with the local image viewer. Confirm:

- six panels appear in the correct order;
- raw and denoised titles match their inputs;
- all pole figures share one MRD upper limit;
- AD/TD-RD direction annotations are readable;
- no pole figure, colorbar, mean label, or P10–P90 range is clipped;
- raw and denoised trend plots use identical x/y limits.

- [ ] **Step 3: Fix only observed layout defects and rerun**

If inspection finds clipping or imbalance, adjust figure position, `TileSpacing`, label offsets, or the colorbar placement in `generate_c_axis_paper_figures.m`, rerun the formal test, and repeat Step 2. Do not change scientific data, ODF parameters, color limits, or summary values during visual iteration.

- [ ] **Step 4: Register the additional figures in README**

Add these entries under `05_texture`:

```markdown
- `05_texture/c_axis_pole_figures_raw.png`
- `05_texture/c_axis_pole_figures_denoised.png`
- `05_texture/c_axis_mean_orientation_raw.png`
- `05_texture/c_axis_mean_orientation_denoised.png`
```

State that raw and denoised are separate processing variants, the pole figures share one MRD scale, and the trend figures share one axis range.

- [ ] **Step 5: Run final verification**

Run:

```powershell
git diff --check
Get-ChildItem 'results/mtex_ebsd_comprehensive/05_texture/c_axis_*.png' |
  Select-Object Name,Length,LastWriteTime
```

Then rerun:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath('C:/Users/22069/Documents/MATLAB/mtex-6.1.1'); startup_mtex; addpath('tools/mtex'); test_c_axis_paper_figures"
```

Expected: no `git diff --check` errors, four non-empty outputs, and `test_c_axis_paper_figures passed`.

- [ ] **Step 6: Commit documentation only**

```powershell
git add -- results/mtex_ebsd_comprehensive/README.md
git commit -m "docs: register c-axis paper figures"
```
