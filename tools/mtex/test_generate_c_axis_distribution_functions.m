function test_generate_c_axis_distribution_functions()
%TEST_GENERATE_C_AXIS_DISTRIBUTION_FUNCTIONS Verify registered artifacts.

textureDir = string(tempname);
inputDir = fullfile(textureDir,"input");
mkdir(textureDir);
mkdir(inputDir);
cleanupObject = onCleanup(@() rmdir(textureDir,"s"));
inputPath = fullfile(inputDir,"synthetic_c_axis_rows.csv");
writetable(build_synthetic_rows(),inputPath);

metadata = generate_c_axis_distribution_functions(textureDir, ...
  struct("input_path",inputPath));
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
  assert(isscalar(info) && info.bytes > 0);
end

ad = readtable(fullfile(textureDir,expected(1)), ...
  "TextType","string","VariableNamingRule","preserve");
spherical = readtable(fullfile(textureDir,expected(2)), ...
  "TextType","string","VariableNamingRule","preserve");
parameters = readtable(fullfile(textureDir,expected(3)), ...
  "TextType","string","VariableNamingRule","preserve");
assert(isequal(string(ad.Properties.VariableNames), [
  "sample","diameter_mm","cold_reduction_percent","variant", ...
  "support","weighting","bin_width_deg","bin_lower_deg", ...
  "bin_upper_deg","bin_center_deg","valid_source_count", ...
  "valid_source_weight","observed_probability","pdf_per_degree", ...
  "random_probability","mrd"]));
assert(isequal(string(spherical.Properties.VariableNames), [
  "sample","diameter_mm","cold_reduction_percent","variant", ...
  "support","weighting","kernel_halfwidth_deg","n_mu","n_phi", ...
  "grid_index","theta_ad_deg","phi_about_ad_deg","cell_weight", ...
  "mrd"]));
assert(isequal(string(parameters.Properties.VariableNames), [
  "scope","sample","variant","support","weighting", ...
  "parameter","value","unit","role","definition"]));

assert(numel(unique(ad.sample(ad.support=="raw_full"))) == 6);
assert(numel(unique(ad.sample(ad.support=="denoised_full"))) == 6);
assert(numel(unique(ad.sample(ad.support=="raw_common"))) == 6);
assert(numel(unique(ad.sample(ad.support=="denoised_raw_common"))) == 6);
assert(~any(contains(ad.support( ...
  ad.weighting=="area_weighted_grain_mean"), ...
  "common")));
assert(all(ad.valid_source_count(ad.support=="raw_full" & ...
  ad.weighting=="pixel_weighted")==3));
assert(all(ad.valid_source_count(ad.support=="denoised_full" & ...
  ad.weighting=="pixel_weighted")==3));
assert(all(ad.valid_source_count(ad.support=="raw_common")==2));
assert(all(ad.valid_source_count( ...
  ad.support=="denoised_raw_common")==2));
assert_state_weight(ad,"7d","raw_full","pixel_weighted",4);
assert_state_weight(ad,"7d","denoised_full","pixel_weighted",5);
assert_state_weight(ad,"7d","raw_common","pixel_weighted",3);
assert_state_weight(ad,"7d","denoised_raw_common", ...
  "pixel_weighted",3);
assert_common_ad_oracle(ad);
adSums = groupsummary(ad, ...
  ["sample","variant","support","weighting","bin_width_deg"], ...
  "sum","observed_probability");
assert(all(abs(adSums.sum_observed_probability-1) < 1e-10));
assert(all(spherical.mrd >= 0 & isfinite(spherical.mrd)));
assert(all(ismember([1;2;5],unique(ad.bin_width_deg))));
assert(all(ismember([5;7.5;10], ...
  unique(spherical.kernel_halfwidth_deg))));
assert(any(parameters.parameter=="peak_grid_distance_deg"));
assert(any(parameters.parameter=="peak_grid_stable"));
requiredParameters = [
  "spherical_kernel_family"
  "theta_ad_definition"
  "phi_about_ad_definition"
  "phi_positive_direction"
  "primary_weighting"
  "primary_supports"
  "common_supports"
  "color_percentile"
  "color_rounding_increment_mrd"
  "color_saturation_limit_fraction"
  "peak_tie_exclusion_deg"
  "raw_minimum_mrd_kernel_5_deg"
  "raw_minimum_mrd_kernel_7_5_deg"
  "raw_minimum_mrd_kernel_10_deg"
  "raw_negative_mass_kernel_5_deg"
  "raw_negative_cell_fraction_kernel_5_deg"
  "maximum_negative_mass_kernel_5_deg"
  "nonnegative_projection_scale_kernel_5_deg"
];
assert(all(ismember(requiredParameters,parameters.parameter)));
tieExclusion = parameters.scope=="global" & ...
  parameters.parameter=="peak_tie_exclusion_deg";
assert(nnz(tieExclusion) == 1 && ...
  double(parameters.value(tieExclusion)) == 10);
negativeMassRows = startsWith(parameters.parameter, ...
  "raw_negative_mass_kernel_");
assert(all(double(parameters.value(negativeMassRows)) >= 0 & ...
  double(parameters.value(negativeMassRows)) <= 0.01));
maximumMassRows = startsWith(parameters.parameter, ...
  "maximum_negative_mass_kernel_");
assert(all(double(parameters.value(maximumMassRows)) == 0.01));
projectionRows = startsWith(parameters.parameter, ...
  "nonnegative_projection_scale_kernel_");
assert(all(double(parameters.value(projectionRows)) > 0));
stabilityRows = parameters.parameter=="peak_grid_stable";
assert(all(ismember(double(parameters.value(stabilityRows)),[0;1])));
assert(nnz(stabilityRows) == 36);
assert_peak_stability(parameters,"7d","raw_full", ...
  "pixel_weighted",1);
assert_peak_stability(parameters,"6.5d","raw_full", ...
  "pixel_weighted",0);
assert(metadata.raw_state_count == 6);
assert(metadata.denoised_state_count == 6);
fprintf("GENERATE_C_AXIS_DISTRIBUTION_FUNCTIONS_TESTS_OK\n");
clear cleanupObject
end

function rows = build_synthetic_rows()
samples = ["7d";"6.5d";"6d";"5.5d";"5d";"4.5d"];
diameters = [7;6.5;6;5.5;5;4.5];
reductions = [0;13.8;26.5;38.3;49;58.7];
rows = table();
for sampleIndex = 1:numel(samples)
  baseTheta = 12 + 8*(sampleIndex-1);
  rawTheta = baseTheta+[0;8;16];
  rawPhi = [0;55;120];
  rawWeight = [1;2;1];
  if sampleIndex == 1
    rawTheta = [0;0;0];
  elseif sampleIndex == 2
    rawTheta = [45;45;45];
    rawPhi = [-90;-90;90];
    rawWeight = [0.5;0.5;1];
  end
  rows = [rows; make_population(samples(sampleIndex), ...
    diameters(sampleIndex),reductions(sampleIndex),"raw", ...
    "pixel_weighted",[1;2;3],rawTheta,rawPhi,rawWeight)]; %#ok<AGROW>
  rows = [rows; make_population(samples(sampleIndex), ...
    diameters(sampleIndex),reductions(sampleIndex),"raw", ...
    "area_weighted_grain_mean",[101;102], ...
    baseTheta+[4;18],[25;145], ...
    [1;3])]; %#ok<AGROW>
  rows = [rows; make_population(samples(sampleIndex), ...
    diameters(sampleIndex),reductions(sampleIndex),"denoised", ...
    "pixel_weighted",[2;3;4],baseTheta+[6;13;19], ...
    [50;115;175],[2;1;2])]; %#ok<AGROW>
  rows = [rows; make_population(samples(sampleIndex), ...
    diameters(sampleIndex),reductions(sampleIndex),"denoised", ...
    "area_weighted_grain_mean",[101;102], ...
    baseTheta+[7;15],[35;155], ...
    [2;1])]; %#ok<AGROW>
end
end

function assert_state_weight(ad,sample,support,weighting,expected)
selected = ad.sample==sample & ad.support==support & ...
  ad.weighting==weighting;
assert(any(selected));
actual = unique(ad.valid_source_weight(selected));
assert(isscalar(actual) && abs(actual-expected) < 1e-12);
end

function assert_common_ad_oracle(ad)
actual = ad(ad.sample=="6d" & ad.variant=="raw" & ...
  ad.support=="raw_common" & ad.weighting=="pixel_weighted" & ...
  ad.bin_width_deg==2,:);
actual = sortrows(actual,"bin_center_deg");
expected = compute_c_axis_ad_distribution( ...
  [36;44],[2;1],0:2:90);
assert(height(actual) == height(expected));
assert(max(abs(actual.observed_probability- ...
  expected.observed_probability)) < 1e-12);
assert(max(abs(actual.pdf_per_degree-expected.pdf_per_degree)) < ...
  1e-12);
assert(max(abs(actual.mrd-expected.mrd)) < 1e-12);
end

function assert_peak_stability( ...
  parameters,sample,support,weighting,expected)
selected = parameters.scope=="state" & ...
  parameters.sample==sample & parameters.support==support & ...
  parameters.weighting==weighting & ...
  parameters.parameter=="peak_grid_stable";
assert(nnz(selected) == 1);
actual = double(parameters.value(selected));
distance = state_parameter_value(parameters,sample,support,weighting, ...
  "peak_grid_distance_deg");
nearPrimary = state_parameter_value( ...
  parameters,sample,support,weighting,"peak_near_tie_primary");
nearAudit = state_parameter_value( ...
  parameters,sample,support,weighting,"peak_near_tie_audit");
assert(actual == expected, ...
  "Peak stability actual=%g expected=%g distance=%g primaryTie=%g auditTie=%g.", ...
  actual,expected,distance,nearPrimary,nearAudit);
end

function value = state_parameter_value( ...
  parameters,sample,support,weighting,parameter)
selected = parameters.scope=="state" & ...
  parameters.sample==sample & parameters.support==support & ...
  parameters.weighting==weighting & ...
  parameters.parameter==parameter;
assert(nnz(selected) == 1);
value = double(parameters.value(selected));
end

function rows = make_population(sample,diameter,reduction,variant, ...
  weighting,sourceId,thetaDeg,phiDeg,sourceWeight)
x = cosd(thetaDeg);
radial = sind(thetaDeg);
y = radial .* cosd(phiDeg);
z = radial .* sind(phiDeg);
n = numel(sourceId);
rows = table(repmat(sample,n,1),repmat(diameter,n,1), ...
  repmat(reduction,n,1),repmat(variant,n,1), ...
  repmat(weighting,n,1),sourceId(:),sourceWeight(:), ...
  sourceWeight(:)/sum(sourceWeight),x(:),y(:),z(:),thetaDeg(:), ...
  phiDeg(:),'VariableNames',cellstr([
  "sample","diameter_mm","cold_reduction_percent","variant", ...
  "weighting","source_id","source_weight","normalized_weight", ...
  "c_axis_x_ad","c_axis_y_td_rd","c_axis_z_nd", ...
  "c_axis_ad_acute_deg","c_axis_azimuth_about_ad_deg"]));
end
