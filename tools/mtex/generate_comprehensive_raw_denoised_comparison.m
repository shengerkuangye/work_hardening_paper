function [comparison, rankAgreement] = ...
  generate_comprehensive_raw_denoised_comparison(outputRoot, options)
%GENERATE_COMPREHENSIVE_RAW_DENOISED_COMPARISON Pair summary metrics.
% Raw values remain primary; denoised values are paired sensitivity results,
% not independent experimental replicates.

arguments
  outputRoot (1,1) string
  options (1,1) struct = struct()
end

contract = comprehensive_ebsd_output_contract();
artifactDir = fullfile(outputRoot, "07_raw_denoised_comparison");
if ~isfolder(artifactDir)
  mkdir(artifactDir);
end
sourceTables = load_source_tables(outputRoot, options);
moduleNames = string(fieldnames(sourceTables));
comparison = empty_comparison_table();

for moduleIndex = 1:numel(moduleNames)
  moduleName = moduleNames(moduleIndex);
  source = sourceTables.(moduleName);
  assert(istable(source) && ~isempty(source), ...
    "Summary table for module %s is empty.", moduleName);
  requiredKeys = ["sample","diameter_mm", ...
    "cold_reduction_percent","variant"];
  assert(all(ismember(requiredKeys, ...
    string(source.Properties.VariableNames))), ...
    "Module %s lacks registered scan keys.", moduleName);
  [dimensionNames, measureNames] = classify_columns(source);
  for measureIndex = 1:numel(measureNames)
    rows = pair_measure(source, moduleName, measureNames(measureIndex), ...
      dimensionNames);
    comparison = append_rows(comparison, rows);
  end
end

comparison = sortrows(comparison, ...
  {'module','metric','cold_reduction_percent'});
[comparison, rankAgreement] = add_agreement_fields(comparison);
comparison = comparison(:, cellstr( ...
  contract.summaryColumns.raw_denoised_metric_comparison));
rankAgreement = rankAgreement(:, cellstr( ...
  contract.summaryColumns.raw_denoised_rank_agreement));
writetable(comparison, fullfile(artifactDir, ...
  "raw_denoised_metric_comparison.csv"));
writetable(rankAgreement, fullfile(artifactDir, ...
  "raw_denoised_rank_agreement.csv"));
if ~isfield(options, "write_figure") || options.write_figure
  plot_rank_agreement(rankAgreement, fullfile(artifactDir, ...
    "raw_denoised_metric_comparison.png"));
end
end

function sourceTables = load_source_tables(outputRoot, options)
if isfield(options, "source_tables")
  sourceTables = options.source_tables;
  assert(isstruct(sourceTables) && isscalar(sourceTables));
  return
end
definitions = { ...
  "grain_morphology", "02_grain_morphology", ...
    "grain_morphology_summary.csv"; ...
  "boundaries", "03_boundaries", "boundary_summary.csv"; ...
  "intragranular", "04_intragranular", ...
    "intragranular_summary.csv"; ...
  "texture", "05_texture", "texture_summary.csv"; ...
  "axial_propensity", "06_axial_propensity", ...
    "axial_propensity_summary.csv"};
sourceTables = struct();
for definitionIndex = 1:size(definitions, 1)
  path = fullfile(outputRoot, definitions{definitionIndex, 2}, ...
    definitions{definitionIndex, 3});
  assert(isfile(path), "Missing summary table: %s", path);
  sourceTables.(definitions{definitionIndex, 1}) = readtable(path, ...
    "TextType", "string", "VariableNamingRule", "preserve");
end
end

function [dimensionNames, measureNames] = classify_columns(source)
names = string(source.Properties.VariableNames);
commonNames = ["sample","diameter_mm","cold_reduction_percent", ...
  "variant"];
parameterNames = ["grain_detection_deg","min_grain_pixels", ...
  "classification_deg","kam_order","kam_threshold_deg", ...
  "texture_kernel_halfwidth_deg","texture_grid_resolution_deg"];
dimensionNames = strings(0,1);
measureNames = strings(0,1);
for name = names
  if ismember(name, commonNames)
    continue
  end
  values = source.(name);
  if isstring(values) || iscellstr(values) || iscategorical(values) || ...
      ismember(name, parameterNames)
    dimensionNames(end+1,1) = name; %#ok<AGROW>
  elseif isnumeric(values) || islogical(values)
    measureNames(end+1,1) = name; %#ok<AGROW>
  end
end
assert(~isempty(measureNames), ...
  "No numeric summary metrics are available for comparison.");
end

function rows = pair_measure(source, moduleName, measureName, dimensions)
metricLabels = metric_label(source, measureName, dimensions);
keys = table(source.sample, source.diameter_mm, ...
  source.cold_reduction_percent, metricLabels, ...
  'VariableNames', {'sample','diameter_mm', ...
  'cold_reduction_percent','metric'});
pairKeys = unique(keys, "rows", "stable");
rows = empty_comparison_table();
for keyIndex = 1:height(pairKeys)
  sameKey = source.sample == pairKeys.sample(keyIndex) & ...
    source.diameter_mm == pairKeys.diameter_mm(keyIndex) & ...
    source.cold_reduction_percent == ...
      pairKeys.cold_reduction_percent(keyIndex) & ...
    metricLabels == pairKeys.metric(keyIndex);
  rawIndex = find(sameKey & string(source.variant) == "raw");
  denoisedIndex = find(sameKey & string(source.variant) == "denoised");
  assert(numel(rawIndex) == 1 && numel(denoisedIndex) == 1, ...
    "Expected one raw/denoised pair for %s %s %s.", ...
    moduleName, measureName, pairKeys.sample(keyIndex));
  rawValue = double(source.(measureName)(rawIndex));
  denoisedValue = double(source.(measureName)(denoisedIndex));
  absoluteDifference = denoisedValue - rawValue;
  if rawValue == 0
    relativeDifference = NaN;
  else
    relativeDifference = 100 * absoluteDifference / abs(rawValue);
  end
  row = table(pairKeys.sample(keyIndex), ...
    pairKeys.diameter_mm(keyIndex), ...
    pairKeys.cold_reduction_percent(keyIndex), moduleName, ...
    pairKeys.metric(keyIndex), rawValue, denoisedValue, ...
    absoluteDifference, relativeDifference, false, false, false, ...
    'VariableNames', rows.Properties.VariableNames);
  rows = [rows; row]; %#ok<AGROW>
end
end

function labels = metric_label(source, measureName, dimensions)
labels = repmat(measureName, height(source), 1);
for dimensionName = dimensions'
  values = string(source.(dimensionName));
  labels = labels + "[" + dimensionName + "=" + values + "]";
end
end

function [comparison, agreement] = add_agreement_fields(comparison)
agreement = table(strings(0,1), strings(0,1), zeros(0,1), ...
  zeros(0,1), zeros(0,1), false(0,1), ...
  'VariableNames', {'module','metric','n_states', ...
  'spearman_rank_agreement','adjacent_direction_agreement_fraction', ...
  'interpretation_change_flag'});
groups = unique(comparison(:, {'module','metric'}), "rows", "stable");
for groupIndex = 1:height(groups)
  mask = comparison.module == groups.module(groupIndex) & ...
    comparison.metric == groups.metric(groupIndex);
  indices = find(mask);
  [~, order] = sort(comparison.cold_reduction_percent(indices));
  indices = indices(order);
  raw = comparison.raw_value(indices);
  denoised = comparison.denoised_value(indices);
  finite = isfinite(raw) & isfinite(denoised);
  rho = spearman_rho(raw(finite), denoised(finite));
  rawDirection = sign(diff(raw));
  denoisedDirection = sign(diff(denoised));
  validTransition = isfinite(raw(1:end-1)) & isfinite(raw(2:end)) & ...
    isfinite(denoised(1:end-1)) & isfinite(denoised(2:end));
  directionAgreement = false(size(rawDirection));
  directionAgreement(validTransition) = ...
    rawDirection(validTransition) == denoisedDirection(validTransition);
  comparison.trend_direction_defined(indices(1)) = false;
  comparison.trend_direction_agrees(indices(1)) = false;
  if numel(indices) > 1
    comparison.trend_direction_defined(indices(2:end)) = validTransition;
    comparison.trend_direction_agrees(indices(2:end)) = ...
      directionAgreement;
  end
  signChanged = sign(raw) ~= sign(denoised) & ...
    isfinite(raw) & isfinite(denoised);
  interpretationChanged = any(~directionAgreement(validTransition)) || ...
    any(signChanged);
  comparison.interpretation_change_flag(indices) = ...
    interpretationChanged;
  if ~any(validTransition)
    directionFraction = NaN;
  else
    directionFraction = mean(directionAgreement(validTransition));
  end
  agreement = [agreement; table(groups.module(groupIndex), ...
    groups.metric(groupIndex), nnz(finite), rho, directionFraction, ...
    interpretationChanged, 'VariableNames', ...
    agreement.Properties.VariableNames)]; %#ok<AGROW>
end
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

function plot_rank_agreement(agreement, outputFile)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [100 100 1500 850]);
cleanupFigure = onCleanup(@() close(figureHandle));
values = agreement.spearman_rank_agreement;
bar(values, "FaceColor", [0.2 0.45 0.75]);
yline(0, "k-");
ylim([-1 1]);
xlabel("Registered summary metric");
ylabel("Raw-denoised Spearman rank agreement");
title(["Raw/denoised paired sensitivity", ...
  " (raw primary; pairs are not replicates)"]);
grid on;
exportgraphics(figureHandle, char(outputFile), "Resolution", 240, ...
  "BackgroundColor", "white");
clear cleanupFigure
end

function output = append_rows(output, rows)
if isempty(output)
  output = rows;
else
  output = [output; rows];
end
end

function result = empty_comparison_table()
result = table(strings(0,1), zeros(0,1), zeros(0,1), ...
  strings(0,1), strings(0,1), zeros(0,1), zeros(0,1), ...
  zeros(0,1), zeros(0,1), false(0,1), false(0,1), false(0,1), ...
  'VariableNames', {'sample','diameter_mm','cold_reduction_percent', ...
  'module','metric','raw_value','denoised_value', ...
  'absolute_difference','relative_difference_percent', ...
  'trend_direction_defined','trend_direction_agrees', ...
  'interpretation_change_flag'});
end
