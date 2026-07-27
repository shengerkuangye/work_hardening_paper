function [sampleSummary,sectionSummary] = ...
  generate_odf_diameter_montage(scanRoot,outputRoot)
%GENERATE_ODF_DIAMETER_MONTAGE Calculate raw alpha-Ti ODF diagnostics.

arguments
  scanRoot (1,1) string
  outputRoot (1,1) string
end

assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded.");

diameterOrder = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"];
phi2Deg = (0:10:60)';
kernelHalfwidthDeg = 5;
gridResolutionDeg = 5;

catalog = comprehensive_ebsd_catalog(scanRoot);
catalog = catalog(catalog.variant == "raw",:);
[isRegistered,catalogOrder] = ismember(diameterOrder,catalog.sample);
assert(all(isRegistered) && numel(unique(catalogOrder)) == 6);
catalog = catalog(catalogOrder,:);
assert(height(catalog) == 6);
assert(isequal(catalog.sample,diameterOrder));

if ~isfolder(outputRoot)
  mkdir(outputRoot);
end

kernel = SO3DeLaValleePoussinKernel( ...
  "halfwidth",kernelHalfwidthDeg * degree);
sampleCount = height(catalog);
sectionCount = numel(phi2Deg);
odfs = cell(sampleCount,1);
validTiHexOrientationCount = zeros(sampleCount,1);
sectionMaximumMrd = zeros(sampleCount,sectionCount);

phi1Deg = 0:gridResolutionDeg:355;
PhiDeg = 0:gridResolutionDeg:90;
[phi1Grid,PhiGrid] = meshgrid(phi1Deg,PhiDeg);

for scanIndex = 1:sampleCount
  [ebsdFull,~] = load_comprehensive_ebsd_scan(catalog(scanIndex,:));
  tiEbsd = ebsdFull("Ti-Hex");
  tiOrientations = tiEbsd.orientations;
  assert(~isempty(tiOrientations));
  assert(string(tiOrientations.CS.pointGroup) == "6/mmm");
  assert(string(tiOrientations.CS.properGroup.pointGroup) == "622");
  assert(string(tiOrientations.SS.pointGroup) == "1");
  [maxPhi1,maxPhi,maxPhi2] = fundamentalRegionEuler( ...
    tiOrientations.CS,tiOrientations.SS);
  assert(isequal([maxPhi1,maxPhi,maxPhi2],[360 90 60] * degree));

  validTiHexOrientationCount(scanIndex) = numel(tiOrientations);
  rbfOdf = calcDensity(tiOrientations,"kernel",kernel, ...
    "weights",ones(numel(tiOrientations),1),"silent");
  rbfOdf = rbfOdf / double(mean(rbfOdf));
  odf = SO3FunHarmonic(rbfOdf,"bandwidth",kernel.bandwidth);
  odf = odf / double(mean(odf));
  odfs{scanIndex} = odf;

  for sectionIndex = 1:sectionCount
    sectionOrientations = orientation.byEuler( ...
      phi1Grid(:) * degree,PhiGrid(:) * degree, ...
      phi2Deg(sectionIndex) * degree, ...
      tiOrientations.CS,tiOrientations.SS);
    sectionMrd = real(eval(odf,sectionOrientations));
    sectionMrd = sectionMrd(:);
    assert(~isempty(sectionMrd) && all(isfinite(sectionMrd)));
    sectionMaximumMrd(scanIndex,sectionIndex) = max(sectionMrd);
  end

  clear ebsdFull tiEbsd tiOrientations rbfOdf odf
end

assert(all(isfinite(sectionMaximumMrd),"all"));
assert(all(sectionMaximumMrd > 0,"all"));
globalMaximumMrd = max(sectionMaximumMrd,[],"all");
assert(isfinite(globalMaximumMrd) && globalMaximumMrd > 0);

sample = catalog.sample;
diameter_mm = catalog.diameter_mm;
cold_reduction_percent = catalog.cold_reduction_percent;
input_path = catalog.input_path;
valid_ti_hex_orientation_count = validTiHexOrientationCount;
crystal_symmetry = repmat("6/mmm",sampleCount,1);
rotational_fundamental_zone = repmat("622",sampleCount,1);
specimen_symmetry = repmat("1",sampleCount,1);
kernel_halfwidth_deg = repmat(kernelHalfwidthDeg,sampleCount,1);
grid_resolution_deg = repmat(gridResolutionDeg,sampleCount,1);
maximum_section_mrd = max(sectionMaximumMrd,[],2);
global_color_limit_max_mrd = repmat(globalMaximumMrd,sampleCount,1);
sampleSummary = table(sample,diameter_mm,cold_reduction_percent, ...
  input_path,valid_ti_hex_orientation_count,crystal_symmetry, ...
  rotational_fundamental_zone,specimen_symmetry, ...
  kernel_halfwidth_deg,grid_resolution_deg,maximum_section_mrd, ...
  global_color_limit_max_mrd);

sample = repelem(catalog.sample,sectionCount);
diameter_mm = repelem(catalog.diameter_mm,sectionCount);
cold_reduction_percent = repelem(cold_reduction_percent,sectionCount);
input_path = repelem(catalog.input_path,sectionCount);
valid_ti_hex_orientation_count = repelem( ...
  validTiHexOrientationCount,sectionCount);
crystal_symmetry = repmat("6/mmm",sampleCount * sectionCount,1);
rotational_fundamental_zone = repmat( ...
  "622",sampleCount * sectionCount,1);
specimen_symmetry = repmat("1",sampleCount * sectionCount,1);
kernel_halfwidth_deg = repmat( ...
  kernelHalfwidthDeg,sampleCount * sectionCount,1);
grid_resolution_deg = repmat( ...
  gridResolutionDeg,sampleCount * sectionCount,1);
phi2_deg = repmat(phi2Deg,sampleCount,1);
section_maximum_mrd = reshape(sectionMaximumMrd.',[],1);
global_color_limit_max_mrd = repmat( ...
  globalMaximumMrd,sampleCount * sectionCount,1);
sectionSummary = table(sample,diameter_mm,cold_reduction_percent, ...
  input_path,valid_ti_hex_orientation_count,crystal_symmetry, ...
  rotational_fundamental_zone,specimen_symmetry, ...
  kernel_halfwidth_deg,grid_resolution_deg,phi2_deg, ...
  section_maximum_mrd,global_color_limit_max_mrd);

assert(height(sampleSummary) == 6);
assert(height(sectionSummary) == 42);
assert(isequal(reshape(sectionSummary.phi2_deg,sectionCount,[]), ...
  repmat(phi2Deg,1,sampleCount)));
assert(globalMaximumMrd == max(sectionSummary.section_maximum_mrd));

writetable(sampleSummary,fullfile(outputRoot,"odf_diameter_summary.csv"));
writetable(sectionSummary, ...
  fullfile(outputRoot,"odf_diameter_section_summary.csv"));
end
