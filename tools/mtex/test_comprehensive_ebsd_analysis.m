function test_comprehensive_ebsd_analysis(projectRoot, outputRoot)
%TEST_COMPREHENSIVE_EBSD_ANALYSIS Verify Task 6 integration invariants.

arguments
  projectRoot (1,1) string = ""
  outputRoot (1,1) string = ""
end

contract = comprehensive_ebsd_output_contract();
test_raw_denoised_pairing(contract);
test_tensile_reliability_and_correlation(contract);
test_manifest_value_serialization();
test_runner_static_safety_and_wording();

if projectRoot ~= ""
  assert(outputRoot ~= "", "outputRoot is required for integration testing.");
  test_full_bundle(string(projectRoot), string(outputRoot), contract);
end

fprintf("test_comprehensive_ebsd_analysis passed\n");
end

function test_manifest_value_serialization()
assert(serialize_manifest_value("raw_full") == "raw_full");
assert(serialize_manifest_value( ...
  ["raw_common","denoised_raw_common"]) == ...
  '["raw_common","denoised_raw_common"]');
assert(serialize_manifest_value([1 2 5]) == "[1 2 5]");
assert(serialize_manifest_value(struct("enabled",true)) == ...
  '{"enabled":true}');
end

function test_raw_denoised_pairing(contract)
samples = repelem(["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"], 2);
diameters = repelem([7;6.48;6.02;5.6;5.25;5], 2);
reductions = repelem([0;14.31;26.04;36;43.75;48.98], 2);
variants = repmat(["raw";"denoised"], 6, 1);
rawMetric = [0;1;1;3;NaN;5];
denoisedMetric = [2;1;2;2;4;NaN];
metric_a = zeros(12,1);
metric_a(1:2:end) = rawMetric;
metric_a(2:2:end) = denoisedMetric;
summary = table(samples, diameters, reductions, variants, metric_a, ...
  'VariableNames', {'sample','diameter_mm','cold_reduction_percent', ...
  'variant','metric_a'});
sourceTables = struct("synthetic", summary);
testRoot = string(tempname);
cleanup = onCleanup(@() remove_test_output(testRoot));
[comparison, agreement] = ...
  generate_comprehensive_raw_denoised_comparison(testRoot, ...
  struct("source_tables", sourceTables, "write_figure", false));

assert(isequal(string(comparison.Properties.VariableNames), ...
  contract.summaryColumns.raw_denoised_metric_comparison));
assert(isequal(string(agreement.Properties.VariableNames), ...
  contract.summaryColumns.raw_denoised_rank_agreement));
assert(height(comparison) == 6);
assert(height(unique(comparison(:, ...
  {'sample','diameter_mm','cold_reduction_percent','module','metric'}), ...
  'rows')) == height(comparison));
finiteDifference = isfinite(comparison.raw_value) & ...
  isfinite(comparison.denoised_value);
assert(all(comparison.denoised_value(finiteDifference) - ...
  comparison.raw_value(finiteDifference) == ...
  comparison.absolute_difference(finiteDifference)));
zeroRaw = comparison.sample == "7d";
assert(nnz(zeroRaw) == 1);
assert(isnan(comparison.relative_difference_percent(zeroRaw)));
assert(comparison.absolute_difference(zeroRaw) == 2);
assert(isequal(comparison.trend_direction_defined, ...
  [false;true;true;true;false;false]));
assert(~any(comparison.trend_direction_agrees));
assert(height(agreement) == 1);
assert(agreement.n_states == 4);
assert(abs(agreement.spearman_rank_agreement) < 1e-12);
assert(agreement.adjacent_direction_agreement_fraction == 0);
assert(agreement.interpretation_change_flag);
assert(isfile(fullfile(testRoot, "07_raw_denoised_comparison", ...
  "raw_denoised_metric_comparison.csv")));
assert(isfile(fullfile(testRoot, "07_raw_denoised_comparison", ...
  "raw_denoised_rank_agreement.csv")));
end

function test_tensile_reliability_and_correlation(contract)
registeredSamples = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"];
diameters = [7;6.48;6.02;5.6;5.25;5];
reductions = [0;14.31;26.04;36;43.75;48.98];
[ebsd, expectedCorrelationKeys, primaryMetric] = ...
  synthetic_registered_ebsd_comparison(contract, registeredSamples, ...
  diameters, reductions);

detail = synthetic_tensile_detail(diameters);
accepted = synthetic_user_summary(diameters);
testRoot = string(tempname);
cleanup = onCleanup(@() remove_test_output(testRoot));
[merged, correlations, aggregates] = ...
  generate_comprehensive_tensile_integration("", testRoot, ebsd, ...
  struct("tensile_detail", detail, "user_summary", accepted, ...
  "write_figure", false));

assert(isequal(string(merged.Properties.VariableNames), ...
  contract.summaryColumns.ebsd_tensile_merged));
assert(isequal(string(correlations.Properties.VariableNames), ...
  contract.summaryColumns.ebsd_tensile_rank_correlations));
assert(isequal(string(aggregates.Properties.VariableNames), ...
  contract.summaryColumns.tensile_condition_aggregates));
assert(all(ismember(unique(detail.yield_status), ...
  unique(merged.yield_status))));
for detailIndex = 1:height(detail)
  rows = merged.tensile_sample == detail.sample(detailIndex);
  assert(any(rows));
  assert(all(merged.yield_status(rows) == ...
    detail.yield_status(detailIndex)));
end
badStatusRp = merged.tensile_sample == "6.48 Y-2" & ...
  merged.tensile_metric == "Rp0.2_MPa";
assert(all(merged.yield_status(badStatusRp) == ...
  "not_reliable_no_elastic_segment_or_offset_crossing"));
assert(~any(merged.included_in_aggregate(badStatusRp)));
assert(all(startsWith(merged.reliability_rule(badStatusRp), "excluded:")));
assert(all(contains(merged.reliability_rule(badStatusRp), ...
  "not_reliable_no_elastic_segment_or_offset_crossing")));
for absentReason = ["not listed in user-accepted repeats","non-finite"]
  assert(~any(contains(merged.reliability_rule(badStatusRp), absentReason)));
end
badStatusOtherMetric = merged.tensile_sample == "6.48 Y-2" & ...
  merged.tensile_metric == "UTS_engineering_MPa";
assert(all(merged.included_in_aggregate(badStatusOtherMetric)), ...
  "A yield-status failure must not discard finite non-yield quantities.");

userRejected = merged.tensile_sample == "6.02 Y-2";
assert(~any(merged.included_in_aggregate(userRejected)));
assert(all(contains(merged.reliability_rule(userRejected), ...
  "not listed in user-accepted repeats")));
userRejectedRp = userRejected & merged.tensile_metric == "Rp0.2_MPa";
for absentReason = ["non-finite","Rp0.2 yield_status"]
  assert(~any(contains(merged.reliability_rule(userRejectedRp), ...
    absentReason)));
end

nonfiniteRp = merged.tensile_sample == "5.6 Y-2" & ...
  merged.tensile_metric == "Rp0.2_MPa";
assert(~any(merged.included_in_aggregate(nonfiniteRp)));
assert(all(contains(merged.reliability_rule(nonfiniteRp), "non-finite")));
for absentReason = ["not listed in user-accepted repeats", ...
    "Rp0.2 yield_status"]
  assert(~any(contains(merged.reliability_rule(nonfiniteRp), absentReason)));
end
includedRows = merged.included_in_aggregate;
assert(all(startsWith(merged.reliability_rule(includedRows), "included:")));
rp602 = aggregates.diameter_mm == 6.02 & ...
  aggregates.tensile_metric == "Rp0.2_MPa";
assert(nnz(rp602) == 1);
assert(aggregates.tensile_n_valid(rp602) == 1);
assert(aggregates.tensile_mean(rp602) == 886);
assert(aggregates.tensile_sd(rp602) == 0);
mergedKey = merged(:, {'sample','ebsd_metric','tensile_sample', ...
  'tensile_repeat','tensile_metric'});
assert(height(unique(mergedKey, 'rows')) == height(merged), ...
  "The registered Cartesian long-table key must be unique.");
assert(all(correlations.n_states == 6));
assert(all(correlations.tensile_n_valid_min == 1));
registeredTensileMetrics = ["Rp0.2_MPa";"UTS_engineering_MPa"; ...
  "strain_at_UTS_uniform_elongation_percent"];
assert(isequal(sort(unique(aggregates.tensile_metric)), ...
  sort(registeredTensileMetrics)), ...
  "Post-necking apparent true quantities must not enter paper-facing summaries.");
assert(isequal(sort(unique(correlations.tensile_metric)), ...
  sort(registeredTensileMetrics)), ...
  "Post-necking apparent true quantities must not enter EBSD correlations.");
for forbiddenMetric = ["valid_pixel_count","grain_count", ...
    "total_area_um2","boundary_touching_grain_fraction", ...
    "total_boundary_length_um"]
  assert(~any(contains(correlations.ebsd_metric, forbiddenMetric)));
end
actualCorrelationKeys = sort(unique(correlations.ebsd_metric));
assert(isequal(actualCorrelationKeys, sort(expectedCorrelationKeys)), ...
  "Correlation output must match the complete registered EBSD whitelist.");
assert(height(correlations) == 3 * numel(expectedCorrelationKeys));
utsCorrelation = correlations.ebsd_metric == ...
  "intragranular." + primaryMetric & ...
  correlations.tensile_metric == "UTS_engineering_MPa";
assert(nnz(utsCorrelation) == 1);
assert(abs(correlations.spearman_rho(utsCorrelation) - 33/35) < 1e-12);
assert(abs(correlations.leave_one_out_rho_min(utsCorrelation) - 0.9) ...
  < 1e-12);
assert(abs(correlations.leave_one_out_rho_max(utsCorrelation) - 1) ...
  < 1e-12);
assert(isfile(fullfile(testRoot, "08_tensile_integration", ...
  "tensile_condition_aggregates.csv")));
assert(all(contains(correlations.interpretation, ...
  "descriptive Spearman")));
assert(~contains(fileread(which( ...
  "generate_comprehensive_tensile_integration")), "fitlm("));
end

function [ebsd, expectedKeys, primaryMetric] = ...
  synthetic_registered_ebsd_comparison(contract, samples, diameters, ...
  reductions)
morphologySuffix = ...
  "[grain_detection_deg=2][min_grain_pixels=5]";
boundarySuffix = ...
  "[grain_detection_deg=2][classification_deg=15]";
intragranularSuffix = "[grain_detection_deg=2]" + ...
  "[min_grain_pixels=5][kam_order=1][kam_threshold_deg=5]";
textureSuffix = "[weighting=pixel_weighted]" + ...
  "[texture_kernel_halfwidth_deg=5][texture_grid_resolution_deg=5]";
families = ["basal_a";"prismatic_a";"pyramidal_a";"pyramidal_ca"; ...
  "extension_twin_t1";"contraction_twin_c1"];
schmidAssumption = repmat( ...
  "geometric absolute Schmid factor; no CRSS activity claim; " + ...
  "Taylor sensitivity not computed because CRSS was not supplied", 6, 1);
isTwinFamily = contains(families, "twin");
schmidAssumption(isTwinFamily) = schmidAssumption(isTwinFamily) + ...
  "; absolute-value twin screen does not resolve polarity";

registeredModules = [ ...
  repmat("grain_morphology",4,1)
  repmat("boundaries",3,1)
  repmat("intragranular",3,1)
  repmat("texture",5,1)
  repmat("axial_propensity",6,1)];
registeredMetrics = [ ...
  "ecd_area_weighted_median_um" + morphologySuffix
  "aspect_ratio_area_weighted_median" + morphologySuffix
  "shape_factor_number_median" + morphologySuffix
  "long_axis_ad_angle_area_weighted_median_deg" + morphologySuffix
  "boundary_line_density_per_um" + boundarySuffix
  "lagb_2_15_length_fraction" + boundarySuffix
  "twin_candidate_length_fraction" + boundarySuffix
  "kam_mean_deg" + intragranularSuffix
  "grod_mean_deg" + intragranularSuffix
  "gos_area_weighted_mean_deg" + intragranularSuffix
  "texture_index" + textureSuffix
  "m_index" + textureSuffix
  "c_axis_ad_mean_deg" + textureSuffix
  "c_axis_within_30deg_fraction" + textureSuffix
  "c_axis_azimuth_resultant" + textureSuffix
  "area_weighted_mean_max_schmid[family=" + families + ...
    "][assumption=" + schmidAssumption + "]"];
assert(numel(registeredMetrics) == 21);

supportModules = [ ...
  repmat("grain_morphology",3,1)
  repmat("boundaries",2,1)
  repmat("intragranular",2,1)
  "texture"
  "axial_propensity"];
supportMetrics = [ ...
  "grain_count" + morphologySuffix
  "total_area_um2" + morphologySuffix
  "boundary_touching_grain_fraction" + morphologySuffix
  "total_boundary_length_um" + boundarySuffix
  "hagb_ge15_length_fraction" + boundarySuffix
  "valid_pixel_count" + intragranularSuffix
  "grod_mean_deg[grain_detection_deg=2]" + ...
    "[min_grain_pixels=5][kam_order=2][kam_threshold_deg=2]"
  "texture_index[weighting=area_weighted_grain_mean]" + ...
    "[texture_kernel_halfwidth_deg=5][texture_grid_resolution_deg=5]"
  "area_weighted_median_max_schmid[family=basal_a]" + ...
    "[assumption=" + schmidAssumption(1) + "]"];

modules = [registeredModules;supportModules];
metrics = [registeredMetrics;supportMetrics];
nDefinitions = numel(metrics);
nConditions = numel(samples);
sample = repmat(samples, nDefinitions, 1);
diameter = repmat(diameters, nDefinitions, 1);
reduction = repmat(reductions, nDefinitions, 1);
module = repelem(modules, nConditions);
metric = repelem(metrics, nConditions);
rawValue = zeros(nDefinitions * nConditions, 1);
for definitionIndex = 1:nDefinitions
  rows = (definitionIndex - 1) * nConditions + (1:nConditions);
  rawValue(rows) = (1:nConditions)' + 100 * (definitionIndex - 1);
end
denoisedValue = rawValue + 0.01;
ebsd = table(sample, diameter, reduction, module, metric, rawValue, ...
  denoisedValue, 0.01 * ones(size(rawValue)), ones(size(rawValue)), ...
  true(size(rawValue)), true(size(rawValue)), false(size(rawValue)), ...
  'VariableNames', cellstr( ...
  contract.summaryColumns.raw_denoised_metric_comparison));
expectedKeys = registeredModules + "." + registeredMetrics;
primaryMetric = "kam_mean_deg" + intragranularSuffix;
end

function detail = synthetic_tensile_detail(diameters)
sample = strings(12,1);
diameter_mm = repelem(diameters, 2);
area_reduction_percent_vs_7mm = repelem( ...
  [0;14.305306;26.04;36;43.75;48.979592], 2);
repeat = repmat([1;2], 6, 1);
UTS_engineering_MPa = repelem([100;200;300;400;600;500], 2) + ...
  repmat([-1;1], 6, 1);
strain_at_UTS_uniform_elongation_percent = ...
  repelem([28;1.4;1.2;1.7;1.4;1.7], 2);
max_true_stress_MPa = UTS_engineering_MPa + 20;
true_strain_at_max_true_stress_percent = ...
  repelem([32;8;5;6;1.5;5], 2);
Rp0_2_MPa = repelem([580;824;886;872;946;909], 2) + ...
  repmat([-2;2], 6, 1);
yield_status = repmat("ok", 12, 1);
for conditionIndex = 1:6
  sample(2*conditionIndex-1) = string(diameters(conditionIndex)) + " Y-1";
  sample(2*conditionIndex) = string(diameters(conditionIndex)) + " Y-2";
end
sample(1:2) = ["7 M-1";"7 M-2"];
sample(5:6) = ["6.02 Y-1";"6.02 Y-2"];
Rp0_2_MPa(5) = 886;
yield_status(4) = ...
  "not_reliable_no_elastic_segment_or_offset_crossing";
Rp0_2_MPa(8) = NaN;
detail = table(sample, diameter_mm, area_reduction_percent_vs_7mm, ...
  repeat, UTS_engineering_MPa, ...
  strain_at_UTS_uniform_elongation_percent, max_true_stress_MPa, ...
  true_strain_at_max_true_stress_percent, yield_status);
detail.('Rp0.2_MPa') = Rp0_2_MPa;
end

function accepted = synthetic_user_summary(diameters)
cold_reduction_percent_reference = ...
  [0;14.305306;26.04;36;43.75;48.979592];
included_samples = strings(6,1);
excluded_samples = strings(6,1);
for conditionIndex = 1:6
  included_samples(conditionIndex) = ...
    string(diameters(conditionIndex)) + " Y-1;" + ...
    string(diameters(conditionIndex)) + " Y-2";
end
included_samples(1) = "7 M-1;7 M-2";
included_samples(3) = "6.02 Y-1";
excluded_samples(3) = "6.02 Y-2";
valid_repeats_n = [2;2;1;2;2;2];
accepted = table(diameters, cold_reduction_percent_reference, ...
  valid_repeats_n, included_samples, excluded_samples, ...
  'VariableNames', {'diameter_mm','cold_reduction_percent_reference', ...
  'valid_repeats_n','included_samples','excluded_samples'});
end

function test_runner_static_safety_and_wording()
toolDir = fileparts(mfilename("fullpath"));
runnerText = fileread(fullfile(toolDir, ...
  "run_comprehensive_ebsd_analysis.m"));
testText = fileread(fullfile(toolDir, ...
  "test_comprehensive_ebsd_analysis.m"));
fullBundleStart = strfind(testText, "function test_full_bundle");
rejectHelperStart = strfind(testText, ...
  "function assert_rejects_output_root");
assert(~isempty(fullBundleStart) && ~isempty(rejectHelperStart));
fullBundleStart = max(fullBundleStart);
rejectHelperStart = max(rejectHelperStart);
assert(fullBundleStart < rejectHelperStart);
fullBundleText = string(testText(fullBundleStart:rejectHelperStart-1));
assert(count(fullBundleText, 'struct("finalize_only", true)') == 2, ...
  "Existing-bundle integration checks must never enter a fresh full run.");
assert(contains(runnerText, "validate_output_root"));
assert(contains(runnerText, "finalize_only"));
assert(contains(runnerText, ".comprehensive_ebsd_owned"));
assert(contains(runnerText, "assert_owned_output_root(outputRoot)"));
assert(contains(runnerText, "if ~finalizeOnly"));
assert(contains(runnerText, "delimitedTextImportOptions"));
assert(contains(runnerText, "importOptions.DataLines = [2 Inf]"));
assert(contains(runnerText, '["results",".codex_tmp"]'));
assert(contains(runnerText, "sha256_file"));
assert(contains(runnerText, "getMTEXpref"));
assert(~contains(runnerText, "fitlm("));
assert(count(string(runnerText), 'rmdir(outputRoot, "s")') == 1);
assert(~contains(runnerText, "rmdir(projectRoot"));
validatePositions = strfind(runnerText, ...
  char("outputRoot = validate_output_root(projectRoot, outputRoot)"));
deletePosition = strfind(runnerText, char('rmdir(outputRoot, "s")'));
assert(~isempty(validatePositions) && isscalar(deletePosition) && ...
  min(validatePositions) < deletePosition);
assert(contains(runnerText, 'for protectedName = ["data","references"]'));
assert(contains(runnerText, "hashesBefore = hash_catalog_inputs(catalog)"));
assert(contains(runnerText, "hashesAfter = hash_catalog_inputs(catalog)"));
assert(contains(runnerText, "verify_manifest_hashes"));
assert(contains(runnerText, '"input_tensile"'));
assert(contains(runnerText, ...
  '"gr4b23271_cold_deformation_tensile_summary.csv"'));
assert(contains(runnerText, ...
  '"gr4b23271_tensile_final_by_diameter_user_exclude.csv"'));
assert(contains(runnerText, "tensileHashesBefore"));
assert(contains(runnerText, "tensileHashesAfter"));
assert(contains(runnerText, "tensileOptions.tensile_detail"));
assert(contains(runnerText, "tensileOptions.user_summary"));
assert(contains(runnerText, "ebsd_metric, tensile_sample, tensile_repeat"));
assert(contains(runnerText, "直接测量"));
assert(contains(runnerText, "代理指标"));
assert(contains(runnerText, "推断"));
assert(contains(runnerText, "不能作为因果关系的证明"));
assert(contains(runnerText, "拉伸相关性指标白名单"));
assert(contains(runnerText, "KAM 采用一阶邻域和 5° 排除阈值"));
assert(contains(runnerText, "De la Vallee Poussin"));
assert(contains(runnerText, "反极性二阶矩张量"));
assert(contains(runnerText, "Taylor 因子未计算"));
assert(contains(runnerText, "虽其他拉伸字段仍为有限值"));
assert(contains(runnerText, "max_true_stress_MPa"));
assert(contains(runnerText, "true_strain_at_max_true_stress_percent"));
textureCall = strfind(runnerText, ...
  "generate_comprehensive_intragranular_texture");
distributionCall = strfind(runnerText, ...
  "generate_c_axis_distribution_functions");
assert(isscalar(textureCall) && isscalar(distributionCall));
assert(textureCall < distributionCall, ...
  "C-axis distributions must be generated after Module 05 source tables.");
assert(count(string(runnerText), ...
  'case "05_texture"') == 0, ...
  "Registered texture diagnostics must not be duplicated as extras.");
end

function test_full_bundle(projectRoot, outputRoot, contract)
assert(isfolder(projectRoot));
scanRoot = fullfile(projectRoot, "data", "ebsd_kpl_250221_7_df", ...
  "scans");
catalog = comprehensive_ebsd_catalog(scanRoot);
hashesBefore = strings(height(catalog),1);
for inputIndex = 1:height(catalog)
  hashesBefore(inputIndex) = sha256_file(catalog.input_path(inputIndex));
end
assert_rejects_output_root(projectRoot, projectRoot);
assert_rejects_output_root(projectRoot, fullfile(projectRoot, "results"));
assert_rejects_output_root(projectRoot, fullfile(projectRoot, "data"));
assert_rejects_output_root(projectRoot, fullfile(projectRoot, "references"));
assert_rejects_output_root(projectRoot, string(tempname));
assert_preserves_unowned_existing_output(projectRoot);
sentinelDirectory = string(fileparts(outputRoot));
if ~isfolder(sentinelDirectory)
  mkdir(sentinelDirectory);
end
sentinelPath = string(tempname(sentinelDirectory)) + ".txt";
sentinelId = fopen(char(sentinelPath), "w");
assert(sentinelId >= 0);
fprintf(sentinelId, "outside explicit outputRoot\n");
fclose(sentinelId);
cleanupSentinel = onCleanup(@() delete_if_file(sentinelPath));
run_comprehensive_ebsd_analysis(projectRoot, outputRoot, ...
  struct("finalize_only", true));
hashesAfter = strings(height(catalog),1);
for inputIndex = 1:height(catalog)
  hashesAfter(inputIndex) = sha256_file(catalog.input_path(inputIndex));
end
assert(isequal(hashesBefore, hashesAfter));
assert(isfile(sentinelPath), ...
  "Runner deleted a sibling outside the explicit outputRoot.");
verify_manifest_hashes(projectRoot, outputRoot, catalog, hashesBefore, ...
  contract);
firstInventory = bundle_file_inventory(outputRoot);
firstManifest = fileread(fullfile(outputRoot, "analysis_manifest.csv"));
for directory = contract.directories'
  assert(isfolder(fullfile(outputRoot, directory)));
end
assert(isfile(fullfile(outputRoot, "README.md")));
assert(isfile(fullfile(outputRoot, "analysis_manifest.csv")));
verify_readme_inventory(outputRoot, contract, catalog);
merged = readtable(fullfile(outputRoot, "08_tensile_integration", ...
  "ebsd_tensile_merged.csv"), "TextType", "string", ...
  "VariableNamingRule", "preserve");
assert(isequal(string(merged.Properties.VariableNames), ...
  contract.summaryColumns.ebsd_tensile_merged));
assert(nnz(merged.diameter_mm == 6.02 & ...
  merged.tensile_n_valid == 1) > 0);
rankAgreement = readtable(fullfile(outputRoot, ...
  "07_raw_denoised_comparison", "raw_denoised_rank_agreement.csv"), ...
  "TextType", "string", "VariableNamingRule", "preserve");
assert(isequal(string(rankAgreement.Properties.VariableNames), ...
  contract.summaryColumns.raw_denoised_rank_agreement));
correlations = readtable(fullfile(outputRoot, ...
  "08_tensile_integration", "ebsd_tensile_rank_correlations.csv"), ...
  "TextType", "string", "VariableNamingRule", "preserve");
assert(isequal(string(correlations.Properties.VariableNames), ...
  contract.summaryColumns.ebsd_tensile_rank_correlations));
assert(all(correlations.n_states == 6));
aggregates = readtable(fullfile(outputRoot, "08_tensile_integration", ...
  "tensile_condition_aggregates.csv"), "TextType", "string", ...
  "VariableNamingRule", "preserve");
assert(isequal(string(aggregates.Properties.VariableNames), ...
  contract.summaryColumns.tensile_condition_aggregates));
assert(nnz(aggregates.diameter_mm == 6.02 & ...
  aggregates.tensile_n_valid == 1) == 3);
assert(~contains(fileread(fullfile(outputRoot, "README.md")), ...
  "证明了"));
run_comprehensive_ebsd_analysis(projectRoot, outputRoot, ...
  struct("finalize_only", true));
assert(isfile(sentinelPath), ...
  "Repeated run deleted a sibling outside the explicit outputRoot.");
assert(isequal(firstInventory, bundle_file_inventory(outputRoot)), ...
  "Repeated run changed the generated artifact inventory.");
assert(strcmp(firstManifest, fileread(fullfile(outputRoot, ...
  "analysis_manifest.csv"))), ...
  "Repeated run changed deterministic manifest provenance.");
verify_manifest_hashes(projectRoot, outputRoot, catalog, hashesBefore, ...
  contract);
verify_readme_inventory(outputRoot, contract, catalog);
clear cleanupSentinel
end

function assert_rejects_output_root(projectRoot, invalidOutputRoot)
didReject = false;
try
  run_comprehensive_ebsd_analysis(projectRoot, invalidOutputRoot);
catch exception
  didReject = contains(string(exception.message), "outputRoot") || ...
    contains(string(exception.message), "derived-output") || ...
    contains(string(exception.message), "protected source folder");
end
assert(didReject, "Unsafe outputRoot was not rejected: %s", ...
  invalidOutputRoot);
end

function assert_preserves_unowned_existing_output(projectRoot)
resultsRoot = fullfile(projectRoot, "results");
unownedRoot = string(tempname(resultsRoot));
mkdir(unownedRoot);
cleanupUnowned = onCleanup(@() remove_test_output(unownedRoot));
sentinelPath = fullfile(unownedRoot, "must_survive.txt");
fileId = fopen(char(sentinelPath), "w");
assert(fileId >= 0);
fprintf(fileId, "unowned output must survive\n");
fclose(fileId);
didReject = false;
try
  run_comprehensive_ebsd_analysis(projectRoot, unownedRoot);
catch exception
  didReject = contains(string(exception.message), "unowned");
end
assert(didReject && isfile(sentinelPath), ...
  "Runner must reject and preserve an unowned existing outputRoot.");
clear cleanupUnowned
end

function verify_manifest_hashes(projectRoot, outputRoot, catalog, ...
  expectedHashes, contract)
importOptions = delimitedTextImportOptions("NumVariables", 6);
importOptions.DataLines = [2 Inf];
importOptions.Delimiter = ",";
importOptions.VariableNames = ...
  {'category','name','value','path','sha256','evidence_class'};
importOptions.VariableTypes = repmat("string", 1, 6);
importOptions.ExtraColumnsRule = "error";
importOptions.EmptyLineRule = "read";
manifest = readtable(fullfile(outputRoot, "analysis_manifest.csv"), ...
  importOptions);
inputRows = manifest.category == "input_ctf";
assert(nnz(inputRows) == height(catalog));
for inputIndex = 1:height(catalog)
  expectedName = catalog.sample(inputIndex) + "_" + ...
    catalog.variant(inputIndex);
  row = inputRows & manifest.name == expectedName;
  assert(nnz(row) == 1);
  assert(manifest.path(row) == catalog.input_path(inputIndex));
  assert(manifest.sha256(row) == expectedHashes(inputIndex));
end
parameterRows = manifest.category == "registered_parameter";
assert(nnz(parameterRows) == numel(fieldnames(contract.parameters)));
tensileNames = ["tensile_detail";"tensile_user_exclusion_summary"];
tensileRelativePaths = [ ...
  "data/tensile_data/gr4b23271_cold_deformation_2_raw_csv/" + ...
    "gr4b23271_cold_deformation_tensile_summary.csv"; ...
  "data/tensile_data/gr4b23271_cold_deformation_2_raw_csv/" + ...
    "gr4b23271_tensile_final_by_diameter_user_exclude.csv"];
tensileRows = manifest.category == "input_tensile";
assert(nnz(tensileRows) == 2, ...
  "Manifest must contain both tensile integration inputs.");
for tensileIndex = 1:2
  row = tensileRows & manifest.name == tensileNames(tensileIndex);
  assert(nnz(row) == 1);
  assert(manifest.path(row) == tensileRelativePaths(tensileIndex));
  absolutePath = fullfile(projectRoot, ...
    replace(tensileRelativePaths(tensileIndex), "/", filesep));
  assert(manifest.sha256(row) == sha256_file(absolutePath));
end
end

function inventory = bundle_file_inventory(outputRoot)
entries = dir(fullfile(outputRoot, "**", "*"));
entries = entries(~[entries.isdir]);
rootPrefix = string(outputRoot) + string(filesep);
inventory = strings(numel(entries),1);
for entryIndex = 1:numel(entries)
  inventory(entryIndex) = erase(string(fullfile( ...
    entries(entryIndex).folder, entries(entryIndex).name)), rootPrefix);
end
inventory = sort(inventory);
end

function verify_readme_inventory(outputRoot, contract, catalog)
readme = string(fileread(fullfile(outputRoot, "README.md")));
assert(contains(readme, ...
  "(sample, ebsd_metric, tensile_sample, tensile_repeat, tensile_metric)"));
assert(contains(readme, "included_in_aggregate"));
assert(contains(readme, "不能作为独立实验观测或额外重复数"));
for rootArtifact = contract.artifacts.root'
  assert(contains(readme, "`" + rootArtifact + "`"));
end
artifactFields = string(fieldnames(contract.artifacts));
for directoryIndex = 1:numel(contract.directories)
  directory = contract.directories(directoryIndex);
  if directory == "01_standard_maps"
    for catalogIndex = 1:height(catalog)
      path = directory + "/" + catalog.sample(catalogIndex) + "_" + ...
        catalog.variant(catalogIndex) + "_maps.png";
      assert(contains(readme, "`" + path + "`"));
    end
  else
    fieldName = artifactFields(directoryIndex + 1);
    for artifact = string(contract.artifacts.(fieldName))'
      assert(contains(readme, "`" + directory + "/" + artifact + "`"));
    end
  end
end
extraArtifacts = [ ...
  "02_grain_morphology/grain_morphology_quantiles.csv"
  "03_boundaries/boundary_detection_sensitivity.csv"
  "05_texture/texture_numerical_diagnostics.csv"
  "05_texture/raw_denoised_odf_distance.csv"];
for path = extraArtifacts'
  assert(contains(readme, "`" + path + "`"));
end
end

function delete_if_file(path)
if isfile(path)
  delete(path);
end
end

function remove_test_output(outputRoot)
if isfolder(outputRoot)
  rmdir(outputRoot, "s");
end
end
