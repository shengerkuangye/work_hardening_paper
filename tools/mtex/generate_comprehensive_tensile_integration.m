function [merged, correlations, aggregates] = ...
  generate_comprehensive_tensile_integration(projectRoot, outputRoot, ...
  ebsdComparison, options)
%GENERATE_COMPREHENSIVE_TENSILE_INTEGRATION Join raw EBSD and tensile data.
% Correlations are descriptive condition-level associations. They are not
% causal evidence and are not predictive-model fits.

arguments
  projectRoot (1,1) string
  outputRoot (1,1) string
  ebsdComparison table = table()
  options (1,1) struct = struct()
end

contract = comprehensive_ebsd_output_contract();
artifactDir = fullfile(outputRoot, "08_tensile_integration");
if ~isfolder(artifactDir)
  mkdir(artifactDir);
end
if isempty(ebsdComparison)
  comparisonPath = fullfile(outputRoot, ...
    "07_raw_denoised_comparison", ...
    "raw_denoised_metric_comparison.csv");
  assert(isfile(comparisonPath), ...
    "Raw/denoised comparison must be generated before tensile integration.");
  ebsdComparison = readtable(comparisonPath, "TextType", "string", ...
    "VariableNamingRule", "preserve");
end
[detail, userSummary] = load_tensile_sources(projectRoot, options);
detail = normalize_tensile_detail(detail);
userSummary = normalize_user_summary(userSummary);
validate_condition_keys(ebsdComparison, detail, userSummary);

aggregates = aggregate_tensile(detail, userSummary);
aggregates = aggregates(:, cellstr( ...
  contract.summaryColumns.tensile_condition_aggregates));
merged = merge_replicates(ebsdComparison, detail, userSummary, aggregates);
mergedKey = merged(:, {'sample','ebsd_metric','tensile_sample', ...
  'tensile_repeat','tensile_metric'});
assert(height(unique(mergedKey, 'rows')) == height(merged), ...
  "EBSD-tensile Cartesian long-table keys are not unique.");
merged = merged(:, cellstr(contract.summaryColumns.ebsd_tensile_merged));
correlations = correlate_condition_summaries(ebsdComparison, aggregates);
correlations = correlations(:, cellstr( ...
  contract.summaryColumns.ebsd_tensile_rank_correlations));

writetable(merged, fullfile(artifactDir, "ebsd_tensile_merged.csv"));
writetable(correlations, fullfile(artifactDir, ...
  "ebsd_tensile_rank_correlations.csv"));
writetable(aggregates, fullfile(artifactDir, ...
  "tensile_condition_aggregates.csv"));
if ~isfield(options, "write_figure") || options.write_figure
  plot_correlation_summary(correlations, fullfile(artifactDir, ...
    "ebsd_tensile_comparison.png"));
end
end

function [detail, userSummary] = load_tensile_sources(projectRoot, options)
if isfield(options, "tensile_detail")
  detail = options.tensile_detail;
else
  detailPath = fullfile(projectRoot, "data", "tensile_data", ...
    "gr4b23271_cold_deformation_2_raw_csv", ...
    "gr4b23271_cold_deformation_tensile_summary.csv");
  assert(isfile(detailPath), "Missing tensile detail: %s", detailPath);
  detail = readtable(detailPath, "TextType", "string", ...
    "VariableNamingRule", "preserve");
end
if isfield(options, "user_summary")
  userSummary = options.user_summary;
else
  userSummaryPath = fullfile(projectRoot, "data", "tensile_data", ...
    "gr4b23271_cold_deformation_2_raw_csv", ...
    "gr4b23271_tensile_final_by_diameter_user_exclude.csv");
  assert(isfile(userSummaryPath), ...
    "Missing tensile user-exclusion summary: %s", userSummaryPath);
  userSummary = readtable(userSummaryPath, "TextType", "string", ...
    "VariableNamingRule", "preserve");
end
assert(istable(detail) && istable(userSummary));
end

function detail = normalize_tensile_detail(detail)
required = ["sample","diameter_mm","area_reduction_percent_vs_7mm", ...
  "repeat","UTS_engineering_MPa", ...
  "strain_at_UTS_uniform_elongation_percent","Rp0.2_MPa", ...
  "yield_status"];
assert(all(ismember(required, string(detail.Properties.VariableNames))), ...
  "Tensile detail lacks registered fields.");
detail.sample = string(detail.sample);
detail.yield_status = string(detail.yield_status);
assert(all(strlength(detail.yield_status) > 0), ...
  "Every tensile repeat must retain its source yield_status.");
detail.registered_reduction_percent = round( ...
  double(detail.area_reduction_percent_vs_7mm), 2);
end

function userSummary = normalize_user_summary(userSummary)
required = ["diameter_mm","cold_reduction_percent_reference", ...
  "valid_repeats_n","included_samples","excluded_samples"];
assert(all(ismember(required, ...
  string(userSummary.Properties.VariableNames))), ...
  "User-exclusion table lacks registered fields.");
userSummary.included_samples = string(userSummary.included_samples);
userSummary.excluded_samples = string(userSummary.excluded_samples);
userSummary.registered_reduction_percent = round( ...
  double(userSummary.cold_reduction_percent_reference), 2);
assert(height(unique(userSummary(:, ...
  {'diameter_mm','registered_reduction_percent'}), 'rows')) == ...
  height(userSummary), "User-exclusion condition keys are not unique.");
end

function validate_condition_keys(ebsdComparison, detail, userSummary)
requiredEbsd = ["sample","diameter_mm","cold_reduction_percent", ...
  "module","metric","raw_value"];
assert(all(ismember(requiredEbsd, ...
  string(ebsdComparison.Properties.VariableNames))));
ebsdConditions = unique(ebsdComparison(:, ...
  {'sample','diameter_mm','cold_reduction_percent'}), 'rows');
assert(height(ebsdConditions) == 6, ...
  "Expected exactly six EBSD condition keys.");
assert(height(unique(ebsdComparison(:, ...
  {'sample','diameter_mm','cold_reduction_percent','module','metric'}), ...
  'rows')) == height(ebsdComparison), ...
  "Raw EBSD metric keys must be unique.");
for conditionIndex = 1:height(ebsdConditions)
  diameter = ebsdConditions.diameter_mm(conditionIndex);
  reduction = ebsdConditions.cold_reduction_percent(conditionIndex);
  detailRows = detail.diameter_mm == diameter & ...
    detail.registered_reduction_percent == reduction;
  summaryRows = userSummary.diameter_mm == diameter & ...
    userSummary.registered_reduction_percent == reduction;
  assert(any(detailRows) && nnz(summaryRows) == 1, ...
    "Incomplete exact diameter/reduction join for %.4g mm, %.2f%%.", ...
    diameter, reduction);
  includedNames = split_sample_list( ...
    userSummary.included_samples(summaryRows));
  excludedNames = split_sample_list( ...
    userSummary.excluded_samples(summaryRows));
  assert(isempty(intersect(includedNames, excludedNames)), ...
    "Included and excluded tensile sample lists must be disjoint.");
  conditionNames = unique(detail.sample(detailRows));
  registeredNames = unique([includedNames; excludedNames]);
  assert(isequal(sort(conditionNames), sort(registeredNames)), ...
    "Included/excluded lists must exactly cover each tensile condition.");
end
end

function aggregates = aggregate_tensile(detail, userSummary)
metricNames = tensile_metric_names();
sample = strings(0,1);
diameter_mm = zeros(0,1);
cold_reduction_percent = zeros(0,1);
tensile_metric = strings(0,1);
tensile_mean = zeros(0,1);
tensile_sd = zeros(0,1);
tensile_n_valid = zeros(0,1);
aggregates = table(sample, diameter_mm, cold_reduction_percent, ...
  tensile_metric, tensile_mean, tensile_sd, tensile_n_valid);

for conditionIndex = 1:height(userSummary)
  diameter = userSummary.diameter_mm(conditionIndex);
  reduction = userSummary.registered_reduction_percent(conditionIndex);
  conditionRows = find(detail.diameter_mm == diameter & ...
    detail.registered_reduction_percent == reduction);
  acceptedNames = split_sample_list( ...
    userSummary.included_samples(conditionIndex));
  acceptedByUser = ismember(detail.sample(conditionRows), acceptedNames);
  assert(nnz(acceptedByUser) == ...
    userSummary.valid_repeats_n(conditionIndex), ...
    "User-exclusion replicate count mismatch at %.4g mm.", diameter);
  for metricName = metricNames'
    values = double(detail.(metricName)(conditionRows));
    included = acceptedByUser & isfinite(values);
    if metricName == "Rp0.2_MPa"
      included = included & detail.yield_status(conditionRows) == "ok";
    end
    validValues = values(included);
    assert(~isempty(validValues), ...
      "No accepted values for %s at %.4g mm.", metricName, diameter);
    row = table(condition_sample_name(diameter), diameter, reduction, ...
      metricName, mean(validValues), std(validValues, 0), ...
      numel(validValues), 'VariableNames', ...
      aggregates.Properties.VariableNames);
    aggregates = [aggregates; row]; %#ok<AGROW>
  end
end
end

function merged = merge_replicates(ebsdComparison, detail, ...
  userSummary, aggregates)
merged = empty_merged_table();
metricNames = tensile_metric_names();
for ebsdIndex = 1:height(ebsdComparison)
  diameter = ebsdComparison.diameter_mm(ebsdIndex);
  reduction = ebsdComparison.cold_reduction_percent(ebsdIndex);
  summaryIndex = find(userSummary.diameter_mm == diameter & ...
    userSummary.registered_reduction_percent == reduction);
  assert(numel(summaryIndex) == 1);
  acceptedNames = split_sample_list( ...
    userSummary.included_samples(summaryIndex));
  repeatIndices = find(detail.diameter_mm == diameter & ...
    detail.registered_reduction_percent == reduction);
  ebsdMetric = string(ebsdComparison.module(ebsdIndex)) + "." + ...
    string(ebsdComparison.metric(ebsdIndex));
  for repeatIndex = repeatIndices'
    for metricName = metricNames'
      aggregateIndex = find(aggregates.diameter_mm == diameter & ...
        aggregates.cold_reduction_percent == reduction & ...
        aggregates.tensile_metric == metricName);
      assert(numel(aggregateIndex) == 1);
      tensileValue = double(detail.(metricName)(repeatIndex));
      acceptedByUser = ismember(detail.sample(repeatIndex), acceptedNames);
      finiteValue = isfinite(tensileValue);
      yieldReliable = true;
      if metricName == "Rp0.2_MPa"
        yieldReliable = detail.yield_status(repeatIndex) == "ok";
      end
      included = acceptedByUser && finiteValue && yieldReliable;
      rule = replicate_reliability_rule(acceptedByUser, finiteValue, ...
        metricName, detail.yield_status(repeatIndex), included);
      row = table(string(ebsdComparison.sample(ebsdIndex)), ...
        diameter, reduction, ebsdMetric, ...
        double(ebsdComparison.raw_value(ebsdIndex)), ...
        detail.sample(repeatIndex), double(detail.repeat(repeatIndex)), ...
        metricName, tensileValue, detail.yield_status(repeatIndex), ...
        included, aggregates.tensile_mean(aggregateIndex), ...
        aggregates.tensile_sd(aggregateIndex), ...
        aggregates.tensile_n_valid(aggregateIndex), rule, ...
        'VariableNames', merged.Properties.VariableNames);
      merged = [merged; row]; %#ok<AGROW>
    end
  end
end
end

function correlations = correlate_condition_summaries( ...
  ebsdComparison, aggregates)
correlations = empty_correlation_table();
scientificRows = scientific_correlation_metric_mask(ebsdComparison);
assert(any(scientificRows), ...
  "No registered scientific EBSD metrics are available for correlation.");
ebsdComparison = ebsdComparison(scientificRows, :);
ebsdKeys = unique(ebsdComparison(:, {'module','metric'}), ...
  'rows', 'stable');
metricNames = tensile_metric_names();
for ebsdIndex = 1:height(ebsdKeys)
  ebsdRows = ebsdComparison.module == ebsdKeys.module(ebsdIndex) & ...
    ebsdComparison.metric == ebsdKeys.metric(ebsdIndex);
  ebsdState = ebsdComparison(ebsdRows, :);
  ebsdMetric = string(ebsdKeys.module(ebsdIndex)) + "." + ...
    string(ebsdKeys.metric(ebsdIndex));
  for metricName = metricNames'
    tensileState = aggregates(aggregates.tensile_metric == metricName, :);
    [ebsdState, tensileState] = paired_condition_order( ...
      ebsdState, tensileState);
    valid = isfinite(ebsdState.raw_value) & ...
      isfinite(tensileState.tensile_mean);
    if nnz(valid) ~= 6
      % Do not present an incomplete subset as the registered six-state
      % condition-level association (for example, absent Taylor factors).
      continue
    end
    x = ebsdState.raw_value(valid);
    y = tensileState.tensile_mean(valid);
    nValid = tensileState.tensile_n_valid(valid);
    rho = spearman_rho(x, y);
    loo = nan(numel(x), 1);
    for omittedIndex = 1:numel(x)
      keep = true(numel(x), 1);
      keep(omittedIndex) = false;
      loo(omittedIndex) = spearman_rho(x(keep), y(keep));
    end
    finiteLoo = loo(isfinite(loo));
    if isempty(finiteLoo)
      looMinimum = NaN;
      looMaximum = NaN;
    else
      looMinimum = min(finiteLoo);
      looMaximum = max(finiteLoo);
    end
    interpretation = "descriptive Spearman association across " + ...
      "condition-level states with leave-one-state-out sensitivity; " + ...
      "not causal evidence or a predictive fit";
    nValidMinimum = min(nValid);
    nValidMaximum = max(nValid);
    row = table(ebsdMetric, metricName, nnz(valid), rho, ...
      looMinimum, looMaximum, nValidMinimum, nValidMaximum, ...
      interpretation, 'VariableNames', ...
      correlations.Properties.VariableNames);
    correlations = [correlations; row]; %#ok<AGROW>
  end
end
end

function rule = replicate_reliability_rule(acceptedByUser, finiteValue, ...
  metricName, yieldStatus, included)
if included
  rule = "included: user-accepted repeat with finite value";
  if metricName == "Rp0.2_MPa"
    rule = rule + "; yield_status == ok";
  end
  return
end

reasons = strings(0,1);
if ~acceptedByUser
  reasons(end+1,1) = ...
    "tensile sample is not listed in user-accepted repeats"; %#ok<AGROW>
end
if ~finiteValue
  reasons(end+1,1) = "tensile value is non-finite"; %#ok<AGROW>
end
if metricName == "Rp0.2_MPa" && yieldStatus ~= "ok"
  reasons(end+1,1) = "Rp0.2 yield_status is " + yieldStatus; %#ok<AGROW>
end
assert(~isempty(reasons), ...
  "An excluded tensile row must retain at least one explicit reason.");
rule = "excluded: " + join(reasons, "; ");
end

function keep = scientific_correlation_metric_mask(ebsdComparison)
moduleName = string(ebsdComparison.module);
metricLabel = string(ebsdComparison.metric);
morphologySuffix = "[grain_detection_deg=2][min_grain_pixels=5]";
boundarySuffix = "[grain_detection_deg=2][classification_deg=15]";
intragranularSuffix = "[grain_detection_deg=2]" + ...
  "[min_grain_pixels=5][kam_order=1][kam_threshold_deg=5]";
textureSuffix = "[weighting=pixel_weighted]" + ...
  "[texture_kernel_halfwidth_deg=5][texture_grid_resolution_deg=5]";
families = ["basal_a";"prismatic_a";"pyramidal_a"; ...
  "pyramidal_ca";"extension_twin_t1";"contraction_twin_c1"];
allowedKeys = [ ...
  "grain_morphology.ecd_area_weighted_median_um" + morphologySuffix
  "grain_morphology.aspect_ratio_area_weighted_median" + morphologySuffix
  "grain_morphology.shape_factor_number_median" + morphologySuffix
  "grain_morphology.long_axis_ad_angle_area_weighted_median_deg" + morphologySuffix
  "boundaries.boundary_line_density_per_um" + boundarySuffix
  "boundaries.lagb_2_15_length_fraction" + boundarySuffix
  "boundaries.twin_candidate_length_fraction" + boundarySuffix
  "intragranular.kam_mean_deg" + intragranularSuffix
  "intragranular.grod_mean_deg" + intragranularSuffix
  "intragranular.gos_area_weighted_mean_deg" + intragranularSuffix
  "texture.texture_index" + textureSuffix
  "texture.m_index" + textureSuffix
  "texture.c_axis_ad_mean_deg" + textureSuffix
  "texture.c_axis_within_30deg_fraction" + textureSuffix
  "texture.c_axis_azimuth_resultant" + textureSuffix
  "axial_propensity.area_weighted_mean_max_schmid[family=" + families + "]"];
assert(numel(allowedKeys) == 21 && numel(unique(allowedKeys)) == 21);
% Assumption text documents interpretation limits but is not part of the
% physical metric identity.  Normalize only that terminal annotation; all
% numerical definitions and Schmid-family labels remain exact-key fields.
identityLabel = regexprep(metricLabel, "\[assumption=[^\]]*\]$", "");
identityKeys = moduleName + "." + identityLabel;
keep = ismember(identityKeys, allowedKeys);
selectedKeys = unique(identityKeys(keep));
assert(isequal(sort(selectedKeys), sort(allowedKeys)), ...
  "All 21 exact scientific correlation keys must be present.");
selectedDefinitions = unique(moduleName(keep) + "." + metricLabel(keep));
assert(numel(selectedDefinitions) == 21, ...
  "Each scientific metric identity must have one assumption annotation.");
for forbiddenToken = ["valid_pixel_count","grain_count", ...
    "total_area_um2","total_boundary_length_um", ...
    "boundary_touching_grain_fraction"]
  assert(~any(contains(lower(selectedKeys), forbiddenToken)), ...
    "Scientific correlation whitelist admitted support token %s.", ...
    forbiddenToken);
end
end

function [ebsdState, tensileState] = paired_condition_order( ...
  ebsdState, tensileState)
assert(height(ebsdState) == 6 && height(tensileState) == 6, ...
  "Correlations require six condition-level states.");
ebsdState = sortrows(ebsdState, ...
  {'diameter_mm','cold_reduction_percent'});
tensileState = sortrows(tensileState, ...
  {'diameter_mm','cold_reduction_percent'});
assert(isequal(ebsdState.diameter_mm, tensileState.diameter_mm));
assert(isequal(ebsdState.cold_reduction_percent, ...
  tensileState.cold_reduction_percent));
end

function values = split_sample_list(value)
value = strip(string(value));
if ismissing(value) || strlength(value) == 0
  values = strings(0,1);
else
  values = strip(split(value, ";"));
  values = values(values ~= "");
end
end

function names = tensile_metric_names()
names = ["Rp0.2_MPa";"UTS_engineering_MPa"; ...
  "strain_at_UTS_uniform_elongation_percent"];
end

function name = condition_sample_name(diameter)
registered = [7;6.48;6.02;5.6;5.25;5];
names = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"];
index = find(registered == diameter);
assert(numel(index) == 1);
name = names(index);
end

function rho = spearman_rho(x, y)
x = double(x(:));
y = double(y(:));
if numel(x) < 2 || all(x == x(1)) || all(y == y(1))
  rho = NaN;
  return
end
xRank = average_tied_rank(x);
yRank = average_tied_rank(y);
matrix = corrcoef(xRank, yRank);
rho = matrix(1,2);
end

function ranks = average_tied_rank(values)
[sortedValues, order] = sort(values);
sortedRanks = zeros(size(values));
startIndex = 1;
while startIndex <= numel(values)
  endIndex = startIndex;
  while endIndex < numel(values) && ...
      sortedValues(endIndex + 1) == sortedValues(startIndex)
    endIndex = endIndex + 1;
  end
  sortedRanks(startIndex:endIndex) = mean(startIndex:endIndex);
  startIndex = endIndex + 1;
end
ranks = zeros(size(values));
ranks(order) = sortedRanks;
end

function plot_correlation_summary(correlations, outputFile)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [100 100 1550 900]);
cleanupFigure = onCleanup(@() close(figureHandle));
metricNames = unique(correlations.tensile_metric, "stable");
ebsdNames = unique(correlations.ebsd_metric, "stable");
imageValues = nan(numel(ebsdNames), numel(metricNames));
for rowIndex = 1:height(correlations)
  ebsdIndex = find(ebsdNames == correlations.ebsd_metric(rowIndex));
  tensileIndex = find(metricNames == ...
    correlations.tensile_metric(rowIndex));
  imageValues(ebsdIndex, tensileIndex) = correlations.spearman_rho(rowIndex);
end
imagesc(imageValues, [-1 1]);
colormap(turbo(256));
colorbar;
xticks(1:numel(metricNames));
xticklabels(metricNames);
xtickangle(25);
yticks([]);
xlabel("Tensile quantity");
ylabel("Registered raw-EBSD metrics");
title(["Descriptive condition-level Spearman associations", ...
  " (six states; not causal evidence)"]);
exportgraphics(figureHandle, char(outputFile), "Resolution", 240, ...
  "BackgroundColor", "white");
clear cleanupFigure
end

function result = empty_merged_table()
result = table(strings(0,1), zeros(0,1), zeros(0,1), ...
  strings(0,1), zeros(0,1), strings(0,1), zeros(0,1), ...
  strings(0,1), zeros(0,1), strings(0,1), false(0,1), ...
  zeros(0,1), zeros(0,1), zeros(0,1), strings(0,1), ...
  'VariableNames', {'sample','diameter_mm','cold_reduction_percent', ...
  'ebsd_metric','ebsd_value','tensile_sample','tensile_repeat', ...
  'tensile_metric','tensile_value','yield_status', ...
  'included_in_aggregate','tensile_mean','tensile_sd', ...
  'tensile_n_valid','reliability_rule'});
end

function result = empty_correlation_table()
result = table(strings(0,1), strings(0,1), zeros(0,1), ...
  zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
  strings(0,1), 'VariableNames', {'ebsd_metric','tensile_metric', ...
  'n_states','spearman_rho','leave_one_out_rho_min', ...
  'leave_one_out_rho_max','tensile_n_valid_min', ...
  'tensile_n_valid_max','interpretation'});
end
