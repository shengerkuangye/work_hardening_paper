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
