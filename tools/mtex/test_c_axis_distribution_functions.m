function test_c_axis_distribution_functions(textureDir)
%TEST_C_AXIS_DISTRIBUTION_FUNCTIONS Verify c-axis distribution outputs.

test_compute_c_axis_ad_distribution();
test_compute_c_axis_spherical_distribution();
if nargin == 1
  test_registered_generator(string(textureDir));
end
fprintf("C_AXIS_DISTRIBUTION_FUNCTIONS_TESTS_OK\n");
end

function test_registered_generator(textureDir)
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
  assert(isscalar(info) && info.bytes > 0);
end
assert(metadata.raw_state_count == 6);
assert(metadata.denoised_state_count == 6);

contract = comprehensive_ebsd_output_contract();
ad = readtable(fullfile(textureDir,expected(1)), ...
  "TextType","string","VariableNamingRule","preserve");
spherical = readtable(fullfile(textureDir,expected(2)), ...
  "TextType","string","VariableNamingRule","preserve");
parameters = readtable(fullfile(textureDir,expected(3)), ...
  "TextType","string","VariableNamingRule","preserve");
assert(isequal(string(ad.Properties.VariableNames), ...
  contract.summaryColumns.c_axis_ad_distribution_function));
assert(isequal(string(spherical.Properties.VariableNames), ...
  contract.summaryColumns.c_axis_spherical_distribution_function));
assert(isequal(string(parameters.Properties.VariableNames), ...
  contract.summaryColumns.c_axis_distribution_parameters));
assert(numel(unique(ad.sample(ad.support=="raw_full"))) == 6);
assert(numel(unique(ad.sample(ad.support=="denoised_full"))) == 6);
assert(numel(unique(spherical.sample( ...
  spherical.support=="raw_full"))) == 6);
assert(numel(unique(spherical.sample( ...
  spherical.support=="denoised_full"))) == 6);
probabilitySums = groupsummary(ad, ...
  ["sample","variant","support","weighting","bin_width_deg"], ...
  "sum","observed_probability");
assert(all(abs(probabilitySums.sum_observed_probability-1) < 1e-10));
assert(all(isfinite(spherical.mrd) & spherical.mrd >= 0));
spherical.weighted_mrd = spherical.cell_weight .* spherical.mrd;
meanMrd = groupsummary(spherical, ...
  ["sample","variant","support","weighting", ...
  "kernel_halfwidth_deg"],"sum","weighted_mrd");
assert(all(abs(meanMrd.sum_weighted_mrd-1) < 5e-12));
end
