# C-Axis Orientation Distribution Functions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate auditable one-dimensional c-axis-to-AD PDF/MRD curves and complete spherical c-axis MRD maps for all six raw and six denoised EBSD states.

**Architecture:** Add two pure numerical helpers: one computes exact bin-integrated axial PDF/MRD values, and one evaluates an antipodal spherical kernel density on a registered equal-area hemisphere grid. A separate generator reads the existing per-orientation c-axis CSV, builds full- and common-support datasets, writes registered CSV tables, and renders four manuscript figures. The main comprehensive runner and output contract invoke and verify the generator after Module 05 texture data exist.

**Tech Stack:** MATLAB R2025a, MTEX 6.1.1, MATLAB tables and graphics, `vector3d.calcDensity`, `S2DeLaValleePoussinKernel`, CSV/PNG outputs, PowerShell verification.

## Global Constraints

- Work only in `C:\Users\22069\Documents\GitHub\work_hardening_paper\.worktrees\comprehensive-ebsd`.
- Horizontal/specimen x is bar axial direction AD; specimen y is TD/RD; specimen z is ND.
- Treat the alpha-Ti c-axis as antipodal.
- Use `raw_full` as the primary condition result and render `raw_full` and `denoised_full` in separate main figures.
- Use paired `raw_common` and `denoised_raw_common` persistent-source-ID intersections only as processing sensitivity.
- Use pixel/area weighting for main figures; retain area-weighted grain-mean results only in sensitivity tables.
- Use 2 degree one-dimensional bins as primary, with 1 and 5 degree sensitivities.
- Use a 5 degree spherical de la Vallee Poussin halfwidth as primary, with 7.5 and 10 degree sensitivities.
- Use `nMu=36`, `nPhi=144` as the primary equal-area grid and
  `nMu=72`, `nPhi=288` to audit peak-location stability.
- Do not overwrite CTF inputs or existing registered figures.
- Do not treat pixels, grains, denoised rows, or common-mask rows as independent experimental replicates.
- Stage and commit only paths owned by the current task; preserve the pre-existing dirty worktree.

---

## File Structure

**Create**

- `tools/mtex/compute_c_axis_ad_distribution.m` - pure one-dimensional PDF/MRD calculation.
- `tools/mtex/build_c_axis_equal_area_grid.m` - registered equal-area hemisphere evaluation grid.
- `tools/mtex/compute_c_axis_spherical_distribution.m` - antipodal spherical density calculation.
- `tools/mtex/generate_c_axis_distribution_functions.m` - support assembly, CSV writing, and four figures.
- `tools/mtex/test_c_axis_distribution_functions.m` - synthetic, generator, and artifact tests.

**Modify**

- `tools/mtex/comprehensive_ebsd_output_contract.m` - register seven outputs, table columns, and parameters.
- `tools/mtex/test_comprehensive_ebsd_contract.m` - assert the new artifact and schema contract.
- `tools/mtex/run_comprehensive_ebsd_analysis.m` - call the generator after Module 05 is written.
- `tools/mtex/test_comprehensive_ebsd_analysis.m` - verify runner wiring and final inventory.

**Derived outputs, not staged**

- `results/mtex_ebsd_comprehensive/05_texture/c_axis_ad_distribution_function.csv`
- `results/mtex_ebsd_comprehensive/05_texture/c_axis_spherical_distribution_function.csv`
- `results/mtex_ebsd_comprehensive/05_texture/c_axis_distribution_parameters.csv`
- `results/mtex_ebsd_comprehensive/05_texture/c_axis_ad_distribution_raw.png`
- `results/mtex_ebsd_comprehensive/05_texture/c_axis_ad_distribution_denoised.png`
- `results/mtex_ebsd_comprehensive/05_texture/c_axis_spherical_distribution_raw.png`
- `results/mtex_ebsd_comprehensive/05_texture/c_axis_spherical_distribution_denoised.png`

### Task 1: One-dimensional PDF and random-axis MRD

**Files:**

- Create: `tools/mtex/compute_c_axis_ad_distribution.m`
- Create initially, then extend throughout the plan: `tools/mtex/test_c_axis_distribution_functions.m`

**Interfaces:**

- Consumes: `thetaDeg` numeric vector in `[0,90]`, positive `weights` vector, and strictly increasing `binEdgesDeg` from 0 to 90.
- Produces: `[distribution,audit] = compute_c_axis_ad_distribution(thetaDeg,weights,binEdgesDeg)`.
- `distribution` columns: `bin_lower_deg`, `bin_upper_deg`, `bin_center_deg`, `bin_width_deg`, `observed_weight`, `observed_probability`, `pdf_per_degree`, `random_probability`, `mrd`.
- `audit` fields: `valid_source_count`, `valid_source_weight`, `probability_sum`, `random_probability_sum`.

- [ ] **Step 1: Write the failing numerical tests**

Create the test file with:

```matlab
function test_c_axis_distribution_functions(textureDir)
%TEST_C_AXIS_DISTRIBUTION_FUNCTIONS Verify c-axis distribution outputs.

test_one_dimensional_distribution();
if nargin == 1
  test_registered_generator(string(textureDir));
end
fprintf("C_AXIS_DISTRIBUTION_FUNCTIONS_TESTS_OK\n");
end

function test_one_dimensional_distribution()
edges = 0:2:90;
[parallel,audit] = compute_c_axis_ad_distribution( ...
  [0;0],[1;3],edges);
assert(abs(sum(parallel.observed_probability)-1) < 1e-12);
assert(parallel.observed_probability(1) == 1);
assert(all(parallel.observed_probability(2:end) == 0));
assert(abs(audit.valid_source_weight-4) < 1e-12);
assert(abs(audit.random_probability_sum-1) < 1e-12);

normal = compute_c_axis_ad_distribution([90;90],[1;1],edges);
assert(normal.observed_probability(end) == 1);

n = 200000;
mu = ((1:n)'-0.5)/n;
theta = acosd(mu);
random = compute_c_axis_ad_distribution(theta,ones(n,1),edges);
assert(max(abs(random.mrd-1)) < 0.02);

assert_error(@() compute_c_axis_ad_distribution( ...
  [-0.1;20],[1;1],edges));
assert_error(@() compute_c_axis_ad_distribution( ...
  [10;20],[1;0],edges));
assert_error(@() compute_c_axis_ad_distribution( ...
  [10;20],[1;1],0:2:88));
end

function assert_error(functionHandle)
didError = false;
try
  functionHandle();
catch
  didError = true;
end
assert(didError);
end
```

- [ ] **Step 2: Run the test and verify the intended failure**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath('tools/mtex'); test_c_axis_distribution_functions"
```

Expected: MATLAB exits nonzero because
`compute_c_axis_ad_distribution` is undefined.

- [ ] **Step 3: Implement the pure calculation**

Create `compute_c_axis_ad_distribution.m`:

```matlab
function [distribution,audit] = compute_c_axis_ad_distribution( ...
  thetaDeg,weights,binEdgesDeg)
%COMPUTE_C_AXIS_AD_DISTRIBUTION Axial PDF and random-normalized MRD.

thetaDeg = double(thetaDeg(:));
weights = double(weights(:));
binEdgesDeg = double(binEdgesDeg(:)');
assert(~isempty(thetaDeg) && numel(thetaDeg) == numel(weights));
assert(all(isfinite(thetaDeg) & thetaDeg >= 0 & thetaDeg <= 90));
assert(all(isfinite(weights) & weights > 0));
assert(numel(binEdgesDeg) >= 2 && binEdgesDeg(1) == 0 && ...
  binEdgesDeg(end) == 90 && all(diff(binEdgesDeg) > 0));

binCount = numel(binEdgesDeg)-1;
binIndex = discretize(thetaDeg,binEdgesDeg);
assert(all(~isnan(binIndex)));
observedWeight = accumarray(binIndex,weights,[binCount,1],@sum,0);
totalWeight = sum(weights);
observedProbability = observedWeight / totalWeight;
binLowerDeg = binEdgesDeg(1:end-1)';
binUpperDeg = binEdgesDeg(2:end)';
binWidthDeg = binUpperDeg-binLowerDeg;
binCenterDeg = (binLowerDeg+binUpperDeg)/2;
pdfPerDegree = observedProbability ./ binWidthDeg;
randomProbability = cosd(binLowerDeg)-cosd(binUpperDeg);
assert(all(randomProbability > 0));
mrd = observedProbability ./ randomProbability;

distribution = table(binLowerDeg,binUpperDeg,binCenterDeg, ...
  binWidthDeg,observedWeight,observedProbability,pdfPerDegree, ...
  randomProbability,mrd, ...
  "VariableNames",["bin_lower_deg","bin_upper_deg", ...
  "bin_center_deg","bin_width_deg","observed_weight", ...
  "observed_probability","pdf_per_degree", ...
  "random_probability","mrd"]);
audit = struct();
audit.valid_source_count = numel(thetaDeg);
audit.valid_source_weight = totalWeight;
audit.probability_sum = sum(observedProbability);
audit.random_probability_sum = sum(randomProbability);
assert(abs(audit.probability_sum-1) < 1e-12);
assert(abs(audit.random_probability_sum-1) < 1e-12);
end
```

- [ ] **Step 4: Run the focused test**

Run the Step 2 command again.

Expected: MATLAB prints `C_AXIS_DISTRIBUTION_FUNCTIONS_TESTS_OK` and exits
zero.

- [ ] **Step 5: Commit the independently tested helper**

```powershell
git add -- tools/mtex/compute_c_axis_ad_distribution.m tools/mtex/test_c_axis_distribution_functions.m
git commit -m "feat: compute c-axis axial distribution"
```

### Task 2: Equal-area hemisphere grid and antipodal spherical density

**Files:**

- Create: `tools/mtex/build_c_axis_equal_area_grid.m`
- Create: `tools/mtex/compute_c_axis_spherical_distribution.m`
- Modify: `tools/mtex/test_c_axis_distribution_functions.m`

**Interfaces:**

- Produces: `[gridTable,gridVectors] = build_c_axis_equal_area_grid(nMu,nPhi)`.
- Produces: `[distribution,densityFunction,audit] = compute_c_axis_spherical_distribution(cAxisXYZ,weights,nMu,nPhi,kernelHalfwidthDeg)`.
- The distribution table contains `grid_index`, `theta_ad_deg`, `phi_about_ad_deg`, `cell_weight`, and `mrd`.
- `densityFunction` is a normalized antipodal MTEX `S2Fun` for plotting.
- `audit` contains source count/weight, grid mean MRD, maximum location, and canonicalization count.

- [ ] **Step 1: Add failing spherical tests**

Add these calls after `test_one_dimensional_distribution()`:

```matlab
test_equal_area_grid();
test_spherical_distribution();
```

Add:

```matlab
function test_equal_area_grid()
[grid,gridVectors] = build_c_axis_equal_area_grid(18,72);
assert(height(grid) == 18*72);
assert(length(gridVectors) == height(grid));
assert(abs(sum(grid.cell_weight)-1) < 1e-12);
assert(all(grid.theta_ad_deg > 0 & grid.theta_ad_deg < 90));
assert(all(grid.phi_about_ad_deg >= -180 & ...
  grid.phi_about_ad_deg < 180));
end

function test_spherical_distribution()
xyz = [1 0 0;1 0 0;-1 0 0];
weights = [1;2;3];
[parallel,~,audit] = compute_c_axis_spherical_distribution( ...
  xyz,weights,36,144,5);
assert(abs(sum(parallel.cell_weight .* parallel.mrd)-1) < 5e-3);
assert(audit.maximum_theta_ad_deg < 5);
assert(audit.canonicalized_source_count == 1);

[positive,~,~] = compute_c_axis_spherical_distribution( ...
  [0.4 0.8 0.2],[1],18,72,7.5);
[negative,~,~] = compute_c_axis_spherical_distribution( ...
  [-0.4 -0.8 -0.2],[1],18,72,7.5);
assert(max(abs(positive.mrd-negative.mrd)) < 1e-10);

[grid,gridVectors] = build_c_axis_equal_area_grid(36,144);
randomXYZ = [gridVectors.x(:),gridVectors.y(:),gridVectors.z(:)];
random = compute_c_axis_spherical_distribution( ...
  randomXYZ,grid.cell_weight,36,144,10);
assert(max(abs(random.mrd-1)) < 0.15);
end
```

- [ ] **Step 2: Run the test and verify failure**

Run the Task 1 Step 2 MATLAB command.

Expected: failure because `build_c_axis_equal_area_grid` is undefined.

- [ ] **Step 3: Implement the equal-area grid**

Create `build_c_axis_equal_area_grid.m`:

```matlab
function [gridTable,gridVectors] = build_c_axis_equal_area_grid(nMu,nPhi)
%BUILD_C_AXIS_EQUAL_AREA_GRID Equal-area AD-positive hemisphere grid.

assert(isscalar(nMu) && nMu >= 4 && nMu == fix(nMu));
assert(isscalar(nPhi) && nPhi >= 8 && nPhi == fix(nPhi));
mu = ((1:nMu)-0.5)/nMu;
phiDeg = -180 + ((1:nPhi)-0.5) * (360/nPhi);
[phiGridDeg,muGrid] = meshgrid(phiDeg,mu);
thetaGridDeg = acosd(muGrid);
radial = sqrt(max(0,1-muGrid.^2));
x = muGrid;
y = radial .* cosd(phiGridDeg);
z = radial .* sind(phiGridDeg);
gridVectors = vector3d(x(:),y(:),z(:));
gridIndex = (1:numel(x))';
theta_ad_deg = thetaGridDeg(:);
phi_about_ad_deg = phiGridDeg(:);
cell_weight = repmat(1/numel(x),numel(x),1);
gridTable = table(gridIndex,theta_ad_deg,phi_about_ad_deg, ...
  cell_weight);
end
```

- [ ] **Step 4: Implement the spherical density helper**

Create `compute_c_axis_spherical_distribution.m`:

```matlab
function [distribution,densityFunction,audit] = ...
  compute_c_axis_spherical_distribution( ...
  cAxisXYZ,weights,nMu,nPhi,kernelHalfwidthDeg)
%COMPUTE_C_AXIS_SPHERICAL_DISTRIBUTION Antipodal c-axis MRD.

cAxisXYZ = double(cAxisXYZ);
weights = double(weights(:));
assert(size(cAxisXYZ,2) == 3 && size(cAxisXYZ,1) == numel(weights));
assert(all(isfinite(cAxisXYZ),"all"));
assert(all(isfinite(weights) & weights > 0));
assert(isscalar(kernelHalfwidthDeg) && kernelHalfwidthDeg > 0);
norms = vecnorm(cAxisXYZ,2,2);
assert(all(norms > 0));
cAxisXYZ = cAxisXYZ ./ norms;
tol = 1e-12;
flipMask = cAxisXYZ(:,1) < -tol | ...
  (abs(cAxisXYZ(:,1)) <= tol & cAxisXYZ(:,2) < -tol) | ...
  (abs(cAxisXYZ(:,1)) <= tol & abs(cAxisXYZ(:,2)) <= tol & ...
   cAxisXYZ(:,3) < 0);
cAxisXYZ(flipMask,:) = -cAxisXYZ(flipMask,:);

sourceVectors = vector3d(cAxisXYZ(:,1),cAxisXYZ(:,2),cAxisXYZ(:,3));
[distribution,gridVectors] = build_c_axis_equal_area_grid(nMu,nPhi);
kernel = S2DeLaValleePoussinKernel( ...
  "halfwidth",kernelHalfwidthDeg*degree);
densityFunction = calcDensity(sourceVectors,"weights",weights, ...
  "kernel",kernel,"antipodal");
gridMrd = double(densityFunction.eval(gridVectors));
gridMean = sum(distribution.cell_weight .* gridMrd);
assert(isfinite(gridMean) && gridMean > 0);
densityFunction = densityFunction / gridMean;
gridMrd = gridMrd / gridMean;
distribution.mrd = gridMrd;
[maximumMrd,maximumIndex] = max(gridMrd);

audit = struct();
audit.valid_source_count = size(cAxisXYZ,1);
audit.valid_source_weight = sum(weights);
audit.canonicalized_source_count = nnz(flipMask);
audit.grid_mean_mrd = sum(distribution.cell_weight .* gridMrd);
audit.maximum_mrd = maximumMrd;
audit.maximum_theta_ad_deg = ...
  distribution.theta_ad_deg(maximumIndex);
audit.maximum_phi_about_ad_deg = ...
  distribution.phi_about_ad_deg(maximumIndex);
end
```

- [ ] **Step 5: Run focused tests and MTEX static checks**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath('C:/Users/22069/Documents/MATLAB/mtex-6.1.1'); startup_mtex; addpath('tools/mtex'); test_c_axis_distribution_functions; files={'tools/mtex/build_c_axis_equal_area_grid.m','tools/mtex/compute_c_axis_spherical_distribution.m'}; for k=1:numel(files), assert(isempty(checkcode(files{k},'-id'))); end"
```

Expected: success marker, no assertion failures, MATLAB exit zero.

- [ ] **Step 6: Commit the spherical calculation**

```powershell
git add -- tools/mtex/build_c_axis_equal_area_grid.m tools/mtex/compute_c_axis_spherical_distribution.m tools/mtex/test_c_axis_distribution_functions.m
git commit -m "feat: compute spherical c-axis density"
```

### Task 3: Support assembly, registered CSVs, and four figures

**Files:**

- Create: `tools/mtex/generate_c_axis_distribution_functions.m`
- Modify: `tools/mtex/test_c_axis_distribution_functions.m`

**Interfaces:**

- Consumes: `05_texture/c_axis_orientation_distribution.csv`.
- Produces: `metadata = generate_c_axis_distribution_functions(textureDir,options)`.
- Optional `options.input_path` allows a synthetic generator test.
- Produces the seven exact artifacts registered in the design.

- [ ] **Step 1: Add a failing synthetic generator test**

Add:

```matlab
function test_registered_generator(textureDir)
textureDir = string(textureDir);
assert(isfolder(textureDir));
metadata = generate_c_axis_distribution_functions(textureDir);
expected = [
  "c_axis_ad_distribution_function.csv"
  "c_axis_spherical_distribution_function.csv"
  "c_axis_distribution_parameters.csv"
  "c_axis_ad_distribution_raw.png"
  "c_axis_ad_distribution_denoised.png"
  "c_axis_spherical_distribution_raw.png"
  "c_axis_spherical_distribution_denoised.png"
];
for fileName = expected'
  info = dir(fullfile(textureDir,fileName));
  assert(numel(info) == 1 && info.bytes > 0);
end
ad = readtable(fullfile(textureDir,expected(1)), ...
  "TextType","string");
spherical = readtable(fullfile(textureDir,expected(2)), ...
  "TextType","string");
assert(numel(unique(ad.sample(ad.support=="raw_full"))) == 6);
assert(numel(unique(ad.sample(ad.support=="denoised_full"))) == 6);
assert(all(abs(groupsummary(ad, ...
  ["sample","variant","support","weighting","bin_width_deg"], ...
  "sum","observed_probability").sum_observed_probability-1) < 1e-10));
assert(all(spherical.mrd >= 0 & isfinite(spherical.mrd)));
assert(metadata.raw_state_count == 6);
assert(metadata.denoised_state_count == 6);
end
```

- [ ] **Step 2: Run against a copied temporary Module 05 and verify failure**

Run:

```powershell
$worktree = 'C:\Users\22069\Documents\GitHub\work_hardening_paper\.worktrees\comprehensive-ebsd'
$tempRoot = Join-Path $worktree '.codex_tmp\c_axis_distribution_test'
New-Item -ItemType Directory -Force -Path $tempRoot | Out-Null
Copy-Item -LiteralPath (Join-Path $worktree 'results\mtex_ebsd_comprehensive\05_texture\c_axis_orientation_distribution.csv') -Destination $tempRoot -Force
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath('C:/Users/22069/Documents/MATLAB/mtex-6.1.1'); startup_mtex; addpath('tools/mtex'); test_c_axis_distribution_functions('.codex_tmp/c_axis_distribution_test')"
```

Expected: MATLAB fails because
`generate_c_axis_distribution_functions` is undefined.

- [ ] **Step 3: Implement support construction and table generation**

Create the generator with this public structure:

```matlab
function metadata = generate_c_axis_distribution_functions( ...
  textureDir,options)
%GENERATE_C_AXIS_DISTRIBUTION_FUNCTIONS Write c-axis PDFs and maps.

arguments
  textureDir (1,1) string
  options (1,1) struct = struct()
end
inputPath = fullfile(textureDir,"c_axis_orientation_distribution.csv");
if isfield(options,"input_path")
  inputPath = string(options.input_path);
end
assert(isfile(inputPath));
rows = readtable(inputPath,"TextType","string", ...
  "VariableNamingRule","preserve");
required = ["sample","diameter_mm","cold_reduction_percent", ...
  "variant","weighting","source_id","source_weight", ...
  "c_axis_x_ad","c_axis_y_td_rd","c_axis_z_nd", ...
  "c_axis_ad_acute_deg"];
assert(all(ismember(required,string(rows.Properties.VariableNames))));

supports = build_support_rows(rows);
binWidths = [1 2 5];
kernelHalfwidths = [5 7.5 10];
[adTable,sphericalTable,densityModels,parameterTable,metadata] = ...
  compute_all_outputs(supports,binWidths,kernelHalfwidths);
writetable(adTable,fullfile(textureDir, ...
  "c_axis_ad_distribution_function.csv"));
writetable(sphericalTable,fullfile(textureDir, ...
  "c_axis_spherical_distribution_function.csv"));
writetable(parameterTable,fullfile(textureDir, ...
  "c_axis_distribution_parameters.csv"));
render_ad_figures(adTable,textureDir);
render_spherical_figures(sphericalTable,densityModels, ...
  textureDir,parameterTable);
end
```

Implement `build_support_rows` as a local function that:

1. copies raw rows to `support="raw_full"`;
2. copies denoised rows to `support="denoised_full"`;
3. for each sample and `weighting=="pixel_weighted"`, inner-joins raw and
   denoised rows on `source_id`;
4. writes the paired raw orientations as `raw_common`;
5. writes the paired denoised orientations as `denoised_raw_common`;
6. asserts unique source IDs within each sample/variant/weighting;
7. renormalizes weights within every output population;
8. does not construct grain-mean common supports.

Implement `compute_all_outputs` so every support/weighting state gets:

- one-dimensional rows for bin widths 1, 2, and 5;
- spherical rows for kernel halfwidths 5, 7.5, and 10;
- metadata columns prepended in the exact output-contract order;
- `nMu=36`, `nPhi=144`;
- a primary-grid versus `nMu=72`, `nPhi=288` peak-distance audit at the
  5 degree kernel halfwidth, with a 3 degree stability tolerance;
- a retained primary 5 degree `S2Fun` in `densityModels` for raw-full and
  denoised-full pixel-weighted plotting.

Write the peak-distance result and stability flag to the parameter table. Mark a
spherical maximum in the figure only when the corresponding stability flag is
true; a nearly tied or grid-sensitive maximum remains unmarked.

Construct `parameterTable` with the exact columns:

```matlab
["scope","sample","variant","support","weighting", ...
 "parameter","value","unit","role","definition"]
```

Use `scope="global"` with empty sample/state fields for fixed bins, kernels,
grid sizes, coordinate conventions, and colour rules. Use `scope="state"` for
the primary/audit-grid peak distance, stability flag, colour ceiling, and
saturated-grid fraction.

- [ ] **Step 4: Implement the one-dimensional figures**

Use one local function per figure and this fixed selection:

```matlab
primary = adTable.weighting=="pixel_weighted" & ...
  adTable.bin_width_deg==2;
render_one_ad_figure(adTable(primary & ...
  adTable.support=="raw_full",:),"raw", ...
  fullfile(textureDir,"c_axis_ad_distribution_raw.png"));
render_one_ad_figure(adTable(primary & ...
  adTable.support=="denoised_full",:),"denoised", ...
  fullfile(textureDir,"c_axis_ad_distribution_denoised.png"));
```

`render_one_ad_figure` must:

- sort the six samples by `cold_reduction_percent`;
- use a 2-by-1 `tiledlayout`;
- plot `pdf_per_degree` above and `mrd` below;
- use the same six colours and line styles in both panels;
- plot a neutral horizontal MRD = 1 reference;
- fix x limits to `[0 90]`;
- compute one shared PDF y-limit and one shared MRD y-limit from the union of
  raw-full and denoised-full primary rows before either figure is rendered;
- export at 300 dpi with a white background.

- [ ] **Step 5: Implement the spherical figures**

For the main spherical rows, calculate:

```matlab
pooled = sphericalTable.weighting=="pixel_weighted" & ...
  sphericalTable.kernel_halfwidth_deg==5 & ...
  ismember(sphericalTable.support,["raw_full","denoised_full"]);
colorMaximum = ceil(2*prctile( ...
  sphericalTable.mrd(pooled),99.9))/2;
assert(isfinite(colorMaximum) && colorMaximum > 1);
saturatedFraction = mean(sphericalTable.mrd(pooled) > colorMaximum);
assert(saturatedFraction <= 0.002);
```

Record `colorMaximum` and `saturatedFraction` in the parameter table.

`render_one_spherical_figure` must:

- use a 2-by-3 tiled layout ordered by cold reduction;
- plot each retained primary `S2Fun` with
  `"antipodal","earea","contourf","silent"`;
- use `[0 colorMaximum]` and identical contours for all 12 panels;
- annotate AD, TD/RD, and ND using the existing project convention;
- mark the registered-grid maximum with a small neutral symbol;
- use one shared MRD colour bar;
- export raw and denoised figures separately at 300 dpi.

- [ ] **Step 6: Run the generator test**

Run the Step 2 PowerShell block again.

Expected: the success marker is printed, all seven files have nonzero size, six
raw and six denoised states are present, and MATLAB exits zero.

- [ ] **Step 7: Commit the generator**

```powershell
git add -- tools/mtex/generate_c_axis_distribution_functions.m tools/mtex/test_c_axis_distribution_functions.m
git commit -m "feat: render c-axis distribution functions"
```

### Task 4: Output contract and comprehensive runner integration

**Files:**

- Modify: `tools/mtex/comprehensive_ebsd_output_contract.m`
- Modify: `tools/mtex/test_comprehensive_ebsd_contract.m`
- Modify: `tools/mtex/run_comprehensive_ebsd_analysis.m`
- Modify: `tools/mtex/test_comprehensive_ebsd_analysis.m`

**Interfaces:**

- Consumes the seven Task 3 artifacts.
- Produces a comprehensive fresh run that creates and verifies them before
  manifest and README finalization.

- [ ] **Step 1: Add failing contract assertions**

Extend the expected Module 05 artifact list with exactly:

```matlab
"c_axis_ad_distribution_function.csv"
"c_axis_spherical_distribution_function.csv"
"c_axis_distribution_parameters.csv"
"c_axis_ad_distribution_raw.png"
"c_axis_ad_distribution_denoised.png"
"c_axis_spherical_distribution_raw.png"
"c_axis_spherical_distribution_denoised.png"
```

Add expected schema arrays:

```matlab
expectedSummaryColumns.c_axis_ad_distribution_function = [common, ...
  "support","weighting","bin_width_deg","bin_lower_deg", ...
  "bin_upper_deg","bin_center_deg","valid_source_count", ...
  "valid_source_weight","observed_probability","pdf_per_degree", ...
  "random_probability","mrd"];
expectedSummaryColumns.c_axis_spherical_distribution_function = [common, ...
  "support","weighting","kernel_halfwidth_deg","n_mu","n_phi", ...
  "grid_index","theta_ad_deg","phi_about_ad_deg","cell_weight", ...
  "mrd"];
expectedSummaryColumns.c_axis_distribution_parameters = [ ...
  "scope","sample","variant","support","weighting", ...
  "parameter","value","unit","role","definition"];
```

Assert parameters:

```matlab
assert(isequal(contract.parameters.c_axis_ad_bin_widths_deg,[1 2 5]));
assert(contract.parameters.c_axis_ad_primary_bin_width_deg == 2);
assert(isequal( ...
  contract.parameters.c_axis_spherical_halfwidths_deg,[5 7.5 10]));
assert(contract.parameters.c_axis_spherical_primary_halfwidth_deg == 5);
assert(contract.parameters.c_axis_spherical_n_mu == 36);
assert(contract.parameters.c_axis_spherical_n_phi == 144);
assert(contract.parameters.c_axis_spherical_audit_n_mu == 72);
assert(contract.parameters.c_axis_spherical_audit_n_phi == 288);
assert(contract.parameters.c_axis_spherical_peak_tolerance_deg == 3);
```

- [ ] **Step 2: Run the contract test and verify failure**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath('tools/mtex'); test_comprehensive_ebsd_contract"
```

Expected: assertion failure because the artifacts and parameters are not yet in
the contract.

- [ ] **Step 3: Extend the contract**

Append the seven artifact names to `contract.artifacts.dir_05_texture` without
reordering existing artifacts. Add the exact summary-column arrays from Step 1.
Append:

```matlab
contract.parameters.c_axis_ad_bin_widths_deg = [1 2 5];
contract.parameters.c_axis_ad_primary_bin_width_deg = 2;
contract.parameters.c_axis_spherical_halfwidths_deg = [5 7.5 10];
contract.parameters.c_axis_spherical_primary_halfwidth_deg = 5;
contract.parameters.c_axis_spherical_n_mu = 36;
contract.parameters.c_axis_spherical_n_phi = 144;
contract.parameters.c_axis_spherical_audit_n_mu = 72;
contract.parameters.c_axis_spherical_audit_n_phi = 288;
contract.parameters.c_axis_spherical_peak_tolerance_deg = 3;
contract.parameters.c_axis_distribution_primary_weighting = ...
  "pixel_weighted";
contract.parameters.c_axis_distribution_primary_support = "raw_full";
contract.parameters.c_axis_distribution_common_supports = ...
  ["raw_common","denoised_raw_common"];
contract.parameters.c_axis_distribution_color_percentile = 99.9;
```

- [ ] **Step 4: Wire the generator into the runner**

Immediately after:

```matlab
generate_comprehensive_intragranular_texture(scanRoot, outputRoot);
```

add:

```matlab
generate_c_axis_distribution_functions( ...
  fullfile(outputRoot,"05_texture"));
```

Extend `test_comprehensive_ebsd_analysis.m` to assert:

```matlab
runnerText = string(fileread( ...
  fullfile("tools","mtex","run_comprehensive_ebsd_analysis.m")));
textureCall = strfind(runnerText, ...
  "generate_comprehensive_intragranular_texture");
distributionCall = strfind(runnerText, ...
  "generate_c_axis_distribution_functions");
assert(isscalar(textureCall) && isscalar(distributionCall));
assert(textureCall < distributionCall);
```

Also append the seven files to the integration expected-output inventory.

- [ ] **Step 5: Run contract, focused, and runner-wiring tests**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath('C:/Users/22069/Documents/MATLAB/mtex-6.1.1'); startup_mtex; addpath('tools/mtex'); test_comprehensive_ebsd_contract; test_c_axis_distribution_functions('.codex_tmp/c_axis_distribution_test'); test_comprehensive_ebsd_analysis"
```

Expected: all three success paths exit zero and the new artifacts are present in
the contract inventory.

- [ ] **Step 6: Commit integration files only**

```powershell
git add -- tools/mtex/comprehensive_ebsd_output_contract.m tools/mtex/test_comprehensive_ebsd_contract.m tools/mtex/run_comprehensive_ebsd_analysis.m tools/mtex/test_comprehensive_ebsd_analysis.m
git commit -m "feat: register c-axis distribution outputs"
```

### Task 5: Formal execution, visual QA, and scientific audit

**Files:**

- Verify: all five created MATLAB files.
- Verify: all four modified integration files.
- Generate only: the seven `05_texture` artifacts.

**Interfaces:**

- Produces validated local paper figures and numerical tables.
- Does not stage the multi-gigabyte derived result tree.

- [ ] **Step 1: Run all focused static and numerical tests**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath('C:/Users/22069/Documents/MATLAB/mtex-6.1.1'); startup_mtex; addpath('tools/mtex'); test_c_axis_distribution_functions('.codex_tmp/c_axis_distribution_test'); test_comprehensive_ebsd_contract; files={'tools/mtex/compute_c_axis_ad_distribution.m','tools/mtex/build_c_axis_equal_area_grid.m','tools/mtex/compute_c_axis_spherical_distribution.m','tools/mtex/generate_c_axis_distribution_functions.m','tools/mtex/test_c_axis_distribution_functions.m'}; for k=1:numel(files), issues=checkcode(files{k},'-id'); assert(isempty(issues),files{k}); end; disp('C_AXIS_DISTRIBUTION_STATIC_OK');"
```

Expected: both success markers and MATLAB exit zero.

- [ ] **Step 2: Run the complete scratch-bundle integration**

First run a complete scratch-bundle integration:

```powershell
$worktree = (Resolve-Path 'C:\Users\22069\Documents\GitHub\work_hardening_paper\.worktrees\comprehensive-ebsd').Path
$fullTest = Join-Path $worktree '.codex_tmp\c_axis_distribution_full_bundle'
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath('C:/Users/22069/Documents/MATLAB/mtex-6.1.1'); startup_mtex; addpath('tools/mtex'); run_comprehensive_ebsd_analysis(pwd,fullfile(pwd,'.codex_tmp','c_axis_distribution_full_bundle')); test_comprehensive_ebsd_analysis(pwd,fullfile(pwd,'.codex_tmp','c_axis_distribution_full_bundle'))"
```

Expected: the complete contract passes in the scratch bundle, all seven new
artifacts exist under its `05_texture`, and MATLAB exits zero.

- [ ] **Step 3: Run the generator on the registered current bundle**

Resolve the exact owned output directory, verify it is within the worktree, and
record hashes of the existing input table and existing c-axis figures, then run:

```powershell
$worktree = (Resolve-Path 'C:\Users\22069\Documents\GitHub\work_hardening_paper\.worktrees\comprehensive-ebsd').Path
$outputRoot = (Resolve-Path (Join-Path $worktree 'results\mtex_ebsd_comprehensive')).Path
if (-not $outputRoot.StartsWith($worktree,[System.StringComparison]::OrdinalIgnoreCase)) { throw 'Output root escaped worktree' }
if (-not (Test-Path -LiteralPath (Join-Path $outputRoot '.comprehensive_ebsd_owned'))) { throw 'Owned output marker missing' }
$protected = @(
  (Join-Path $outputRoot '05_texture\c_axis_orientation_distribution.csv'),
  (Join-Path $outputRoot '05_texture\c_axis_pole_figures_raw.png'),
  (Join-Path $outputRoot '05_texture\c_axis_pole_figures_denoised.png'),
  (Join-Path $outputRoot '05_texture\c_axis_mean_orientation_raw.png'),
  (Join-Path $outputRoot '05_texture\c_axis_mean_orientation_denoised.png')
) | Where-Object { Test-Path -LiteralPath $_ }
$before = @{}
foreach ($path in $protected) { $before[$path] = (Get-FileHash -Algorithm SHA256 -LiteralPath $path).Hash }
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath('C:/Users/22069/Documents/MATLAB/mtex-6.1.1'); startup_mtex; addpath('tools/mtex'); generate_c_axis_distribution_functions('results/mtex_ebsd_comprehensive/05_texture'); test_c_axis_distribution_functions('results/mtex_ebsd_comprehensive/05_texture')"
foreach ($path in $protected) {
  $after = (Get-FileHash -Algorithm SHA256 -LiteralPath $path).Hash
  if ($before[$path] -ne $after) { throw "Protected artifact changed: $path" }
}
```

Expected: all seven artifacts are regenerated and the formal test exits zero.

- [ ] **Step 4: Audit numerical conservation and support coverage**

Run:

```powershell
$texture = 'C:\Users\22069\Documents\GitHub\work_hardening_paper\.worktrees\comprehensive-ebsd\results\mtex_ebsd_comprehensive\05_texture'
$ad = Import-Csv -LiteralPath (Join-Path $texture 'c_axis_ad_distribution_function.csv')
$spherical = Import-Csv -LiteralPath (Join-Path $texture 'c_axis_spherical_distribution_function.csv')
$ad | Group-Object sample,variant,support,weighting,bin_width_deg | ForEach-Object {
  $sum = ($_.Group | Measure-Object -Property observed_probability -Sum).Sum
  if ([math]::Abs([double]$sum - 1) -gt 1e-9) { throw "AD probability failed: $($_.Name)" }
}
$spherical | Group-Object sample,variant,support,weighting,kernel_halfwidth_deg | ForEach-Object {
  $mean = ($_.Group | ForEach-Object { [double]$_.cell_weight * [double]$_.mrd } | Measure-Object -Sum).Sum
  if ([math]::Abs([double]$mean - 1) -gt 5e-3) { throw "Spherical MRD mean failed: $($_.Name)" }
}
'C_AXIS_DISTRIBUTION_NUMERICAL_AUDIT_OK'
```

Expected: `C_AXIS_DISTRIBUTION_NUMERICAL_AUDIT_OK`.

- [ ] **Step 5: Inspect all four figures**

Open or render:

```text
results/mtex_ebsd_comprehensive/05_texture/c_axis_ad_distribution_raw.png
results/mtex_ebsd_comprehensive/05_texture/c_axis_ad_distribution_denoised.png
results/mtex_ebsd_comprehensive/05_texture/c_axis_spherical_distribution_raw.png
results/mtex_ebsd_comprehensive/05_texture/c_axis_spherical_distribution_denoised.png
```

Require:

- six states are present and ordered 0 to 48.98% reduction;
- PDF and MRD curves are distinguishable and not clipped;
- raw and denoised one-dimensional axes are identical;
- all spherical panels share the same MRD colour scale;
- AD, TD/RD, and ND are correctly oriented;
- titles, legends, maximum markers, and colour bars do not overlap;
- no smoothing or plotting choice hides multimodality.

If visual QA fails, change only figure layout/style in the generator, rerun the
focused test and Step 2, and inspect again. Do not change scientific
normalization, support definitions, bins, kernels, or colour-limit rule during a
layout iteration.

- [ ] **Step 6: Run final repository checks**

Run:

```powershell
git diff --check
git status --short
```

Expected: no whitespace errors. Only task-owned code/spec/plan changes may be
staged or committed; pre-existing user changes and derived results remain
unstaged.

- [ ] **Step 7: Report scientific outputs**

Report:

- whether the one-dimensional distributions move, broaden, split, or exchange
  weight without a mean shift;
- whether the spherical maxima change mainly in polar angle or azimuth about AD;
- whether conclusions agree between raw-full and denoised-full;
- whether the paired common-support sensitivity changes the interpretation;
- which claims remain descriptive because there is one EBSD field per state.
