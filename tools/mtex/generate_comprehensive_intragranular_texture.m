function [intragranularSummary, textureSummary, odfDistance] = ...
  generate_comprehensive_intragranular_texture(scanRoot, outputRoot, options)
%GENERATE_COMPREHENSIVE_INTRAGRANULAR_TEXTURE Export Task 4 results.
% KAM, GROD, and GOS are orientation-gradient/spread proxies. They are not
% direct measurements of dislocation density.

arguments
  scanRoot (1,1) string
  outputRoot (1,1) string
  options (1,1) struct = struct()
end

assert(isfolder(scanRoot), "EBSD scan folder not found: %s", scanRoot);
assert(~isempty(which("EBSD")), ...
  "MTEX is not loaded in the current MATLAB session.");
contract = comprehensive_ebsd_output_contract();
catalog = comprehensive_ebsd_catalog(scanRoot);
catalogRows = selected_catalog_rows(options, height(catalog));
catalog = catalog(catalogRows, :);

if render_odf_only_requested(options)
  textureDir = ensure_directory(outputRoot, "05_texture");
  odfModels = rebuild_texture_odf_models(catalog, contract.parameters);
  verify_rebuilt_texture_indices(odfModels, catalog, textureDir, contract);
  colorLimits = read_registered_texture_color_limits( ...
    textureDir, height(catalog), contract);
  outputName = recovery_output_name(options);
  plot_texture_montage(odfModels, catalog, colorLimits, "odf", ...
    fullfile(textureDir, outputName), contract.parameters);
  intragranularSummary = table();
  textureSummary = table();
  odfDistance = table();
  return
end

intragranularDir = ensure_directory(outputRoot, "04_intragranular");
textureDir = ensure_directory(outputRoot, "05_texture");
pixelPath = fullfile(intragranularDir, "intragranular_by_pixel.csv");
grainPath = fullfile(intragranularDir, "intragranular_by_grain.csv");
cAxisPath = fullfile(textureDir, "c_axis_orientation_distribution.csv");

parameters = contract.parameters;
reconstructionOptions = struct( ...
  "detection_threshold_deg", parameters.primary_grain_detection_deg, ...
  "min_grain_pixels", parameters.min_grain_pixels);
intraOptions = struct( ...
  "grain_detection_deg", parameters.primary_grain_detection_deg, ...
  "min_grain_pixels", parameters.min_grain_pixels, ...
  "kam_orders", parameters.kam_orders, ...
  "kam_thresholds_deg", parameters.kam_thresholds_deg, ...
  "axis_min_grod_deg", parameters.min_boundary_axis_deg);
textureOptionsBase = struct( ...
  "kernel_halfwidth_deg", parameters.texture_kernel_halfwidth_deg, ...
  "grid_resolution_deg", parameters.texture_grid_resolution_deg, ...
  "c_axis_component_zero_tolerance", ...
    parameters.c_axis_component_zero_tolerance, ...
  "min_grain_pixels", parameters.min_grain_pixels);

intragranularSummary = table();
textureSummary = table();
textureDiagnostics = table();
analysisParameters = table();
odfModels = cell(height(catalog), 2);

for scanIndex = 1:height(catalog)
  catalogRow = catalog(scanIndex, :);
  fprintf("INTRAGRANULAR_TEXTURE sample=%s variant=%s\n", ...
    catalogRow.sample, catalogRow.variant);
  [ebsdFull, ~] = load_comprehensive_ebsd_scan(catalogRow);
  [grains, ebsdFull, ~] = reconstruct_comprehensive_grains( ...
    ebsdFull, reconstructionOptions);
  meta = metadata_struct(catalogRow);

  [scanIntra, scanByPixel, scanByGrain] = ...
    compute_intragranular_metrics(ebsdFull, grains, meta, intraOptions);
  assert(isequal(string(scanIntra.Properties.VariableNames), ...
    contract.summaryColumns.intragranular_summary));
  scanByPixel = add_scan_columns(scanByPixel, catalogRow);
  scanByGrain = add_scan_columns(scanByGrain, catalogRow);
  write_scan_table(scanByPixel, pixelPath, scanIndex == 1);
  write_scan_table(scanByGrain, grainPath, scanIndex == 1);
  intragranularSummary = append_rows(intragranularSummary, scanIntra);

  tiEbsd = ebsdFull("Ti-Hex");
  allTiGrains = grains("Ti-Hex");
  allGrainAreas = abs(double(area(allTiGrains)));
  allGrainAreas = allGrainAreas(:);
  assert(~isempty(allTiGrains) && ...
    all(isfinite(allGrainAreas) & allGrainAreas > 0), ...
    "Ti-Hex grains contain invalid areas for %s %s.", ...
    catalogRow.sample, catalogRow.variant);
  retainedGrain = allTiGrains.numPixel >= parameters.min_grain_pixels;
  totalGrainCount = length(allTiGrains);
  retainedGrainCount = nnz(retainedGrain);
  excludedSmallGrainCount = totalGrainCount - retainedGrainCount;
  excludedSmallGrainAreaFraction = ...
    sum(allGrainAreas(~retainedGrain)) / sum(allGrainAreas);
  tiGrains = allTiGrains(retainedGrain);
  assert(~isempty(tiGrains), ...
    "No Ti-Hex grains pass min_grain_pixels for %s %s.", ...
    catalogRow.sample, catalogRow.variant);
  grainAreas = abs(double(area(tiGrains)));
  grainAreas = grainAreas(:);
  assert(all(isfinite(grainAreas) & grainAreas > 0), ...
    "Ti-Hex grains contain invalid areas for %s %s.", ...
    catalogRow.sample, catalogRow.variant);
  assert(~isempty(tiEbsd) && ~isempty(tiGrains), ...
    "No valid Ti-Hex orientations for %s %s.", ...
    catalogRow.sample, catalogRow.variant);
  textureOptions = textureOptionsBase;
  textureOptions.pixel_ids = double(tiEbsd.id(:));
  textureOptions.grain_ids = double(tiGrains.id(:));
  [scanTexture, scanCAxis, scanModel] = compute_texture_metrics( ...
    tiEbsd.orientations, tiGrains.meanOrientation, grainAreas, ...
    meta, textureOptions);
  scanModel.total_ti_hex_grain_count = totalGrainCount;
  scanModel.retained_ti_hex_grain_count = retainedGrainCount;
  scanModel.excluded_small_grain_count = excludedSmallGrainCount;
  scanModel.excluded_small_grain_area_fraction = ...
    excludedSmallGrainAreaFraction;
  assert(isequal(string(scanTexture.Properties.VariableNames), ...
    contract.summaryColumns.texture_summary));
  assert(isequal(string(scanCAxis.Properties.VariableNames), ...
    contract.summaryColumns.c_axis_orientation_distribution));
  assert(all(isfinite(scanTexture{:, 8:end}), "all"));
  write_scan_table(scanCAxis, cAxisPath, scanIndex == 1);
  textureSummary = append_rows(textureSummary, scanTexture);
  textureDiagnostics = append_rows(textureDiagnostics, ...
    entropy_diagnostic_rows(catalogRow, scanModel));
  analysisParameters = append_rows(analysisParameters, ...
    analysis_parameter_row(catalogRow, scanModel, contract));
  odfModels{scanIndex, 1} = scanModel.pixel_weighted_odf;
  odfModels{scanIndex, 2} = ...
    scanModel.area_weighted_grain_mean_odf;

  clear ebsdFull grains tiEbsd tiGrains allTiGrains scanIntra scanByPixel
  clear scanByGrain scanTexture scanCAxis scanModel
end

expectedIntraRows = height(catalog) * ...
  numel(parameters.kam_orders) * numel(parameters.kam_thresholds_deg);
assert(height(intragranularSummary) == expectedIntraRows);
assert(height(textureSummary) == 2 * height(catalog));
assert(isequal(string(textureDiagnostics.Properties.VariableNames), ...
  contract.summaryColumns.texture_numerical_diagnostics));
assert(isequal(string(analysisParameters.Properties.VariableNames), ...
  contract.summaryColumns.task4_analysis_parameters));
conditionKeys = analysisParameters.sample + "|" + analysisParameters.variant;
assert(height(analysisParameters) == height(catalog) && ...
  numel(unique(conditionKeys)) == height(analysisParameters));
writetable(intragranularSummary, fullfile(intragranularDir, ...
  "intragranular_summary.csv"));
writetable(textureSummary, fullfile(textureDir, "texture_summary.csv"));
writetable(textureDiagnostics, fullfile(textureDir, ...
  "texture_numerical_diagnostics.csv"));
writetable(analysisParameters, fullfile(textureDir, ...
  "task4_analysis_parameters.csv"));

odfDistance = raw_denoised_distances(catalog, odfModels, parameters);
assert(isequal(string(odfDistance.Properties.VariableNames), ...
  contract.summaryColumns.raw_denoised_odf_distance));
writetable(odfDistance, fullfile(textureDir, ...
  "raw_denoised_odf_distance.csv"));
texturePlotColorLimits = texture_plot_color_limits(odfModels, parameters);
assert(isequal(string(texturePlotColorLimits.Properties.VariableNames), ...
  contract.summaryColumns.texture_plot_color_limits));
assert(height(texturePlotColorLimits) == 6 && ...
  all(texturePlotColorLimits.scan_count == height(catalog)));
writetable(texturePlotColorLimits, fullfile(textureDir, ...
  "texture_plot_color_limits.csv"));
plot_intragranular_trends(intragranularSummary, ...
  fullfile(intragranularDir, "intragranular_trends.png"));
plot_texture_trends(textureSummary, ...
  fullfile(textureDir, "texture_trends.png"));
plot_texture_montage(odfModels, catalog, texturePlotColorLimits, ...
  "pole", fullfile(textureDir, "pole_figures.png"), parameters);
plot_texture_montage(odfModels, catalog, texturePlotColorLimits, ...
  "ipf", fullfile(textureDir, "inverse_pole_figures.png"), parameters);
plot_texture_montage(odfModels, catalog, texturePlotColorLimits, ...
  "odf", fullfile(textureDir, "odf_sections.png"), parameters);
end

function result = analysis_parameter_row(catalogRow, model, contract)
parameters = contract.parameters;
sample = string(catalogRow.sample);
diameter_mm = double(catalogRow.diameter_mm);
cold_reduction_percent = double(catalogRow.cold_reduction_percent);
variant = string(catalogRow.variant);
grain_detection_deg = double(parameters.primary_grain_detection_deg);
min_grain_pixels = double(parameters.min_grain_pixels);
kam_orders = join(string(parameters.kam_orders), "|");
kam_thresholds_deg = join(string(parameters.kam_thresholds_deg), "|");
grod_axis_cutoff_deg = double(parameters.min_boundary_axis_deg);
grod_axis_aggregation = string(parameters.grod_axis_aggregation);
texture_kernel_family = string(parameters.texture_kernel_family);
texture_kernel_halfwidth_deg = double(parameters.texture_kernel_halfwidth_deg);
harmonic_bandwidth = double(model.harmonic_bandwidth);
texture_evaluation_grid_resolution_deg = ...
  double(parameters.texture_grid_resolution_deg);
coordinate_x = string(parameters.coordinate_x);
coordinate_y = string(parameters.coordinate_y);
coordinate_z = string(parameters.coordinate_z);
primary_variant = string(parameters.primary_variant);
comparison_variant = string(parameters.comparison_variant);
c_axis_hemisphere_convention = ...
  string(parameters.c_axis_hemisphere_convention);
c_axis_component_zero_tolerance = ...
  double(parameters.c_axis_component_zero_tolerance);
clipped_grid_entropy_definition = ...
  string(parameters.clipped_grid_entropy_definition);
odf_harmonic_l2_distance_definition = ...
  string(parameters.raw_denoised_odf_distance_metric);
result = table(sample, diameter_mm, cold_reduction_percent, variant, ...
  grain_detection_deg, min_grain_pixels, kam_orders, kam_thresholds_deg, ...
  grod_axis_cutoff_deg, grod_axis_aggregation, texture_kernel_family, ...
  texture_kernel_halfwidth_deg, harmonic_bandwidth, ...
  texture_evaluation_grid_resolution_deg, coordinate_x, coordinate_y, ...
  coordinate_z, primary_variant, comparison_variant, ...
  c_axis_hemisphere_convention, c_axis_component_zero_tolerance, ...
  clipped_grid_entropy_definition, odf_harmonic_l2_distance_definition);
end

function result = entropy_diagnostic_rows(catalogRow, model)
diagnostics = {model.pixel_weighted_entropy_diagnostics; ...
  model.area_weighted_grain_mean_entropy_diagnostics};
weighting = ["pixel_weighted"; "area_weighted_grain_mean"];
sample = repmat(catalogRow.sample, 2, 1);
diameter_mm = repmat(catalogRow.diameter_mm, 2, 1);
cold_reduction_percent = repmat( ...
  catalogRow.cold_reduction_percent, 2, 1);
variant = repmat(catalogRow.variant, 2, 1);
texture_kernel_halfwidth_deg = repmat(model.kernel_halfwidth_deg, 2, 1);
texture_grid_resolution_deg = repmat(model.grid_resolution_deg, 2, 1);
harmonic_bandwidth = repmat(model.harmonic_bandwidth, 2, 1);
total_ti_hex_grain_count = repmat( ...
  model.total_ti_hex_grain_count, 2, 1);
retained_ti_hex_grain_count = repmat( ...
  model.retained_ti_hex_grain_count, 2, 1);
excluded_small_grain_count = repmat( ...
  model.excluded_small_grain_count, 2, 1);
excluded_small_grain_area_fraction = repmat( ...
  model.excluded_small_grain_area_fraction, 2, 1);
density_normalization_factor = [ ...
  model.pixel_weighted_density_normalization_factor; ...
  model.area_weighted_grain_mean_density_normalization_factor];
grid_point_count = cellfun(@(x) x.grid_point_count, diagnostics);
minimum_raw_odf = cellfun(@(x) x.minimum_raw_odf, diagnostics);
negative_grid_fraction = cellfun( ...
  @(x) x.negative_grid_fraction, diagnostics);
clipped_negative_l1_fraction = cellfun( ...
  @(x) x.clipped_negative_l1_fraction, diagnostics);
clipped_normalization_mean = cellfun( ...
  @(x) x.clipped_normalization_mean, diagnostics);
result = table(sample, diameter_mm, cold_reduction_percent, variant, ...
  weighting, texture_kernel_halfwidth_deg, ...
  texture_grid_resolution_deg, harmonic_bandwidth, ...
  total_ti_hex_grain_count, retained_ti_hex_grain_count, ...
  excluded_small_grain_count, excluded_small_grain_area_fraction, ...
  grid_point_count, ...
  density_normalization_factor, minimum_raw_odf, negative_grid_fraction, ...
  clipped_negative_l1_fraction, clipped_normalization_mean);
end

function requested = render_odf_only_requested(options)
requested = false;
if ~isfield(options, "render_odf_only")
  return
end
value = options.render_odf_only;
assert(isscalar(value) && (islogical(value) || isnumeric(value)) && ...
  isfinite(double(value)) && any(double(value) == [0 1]), ...
  "render_odf_only must be a scalar logical value.");
requested = logical(value);
end

function models = rebuild_texture_odf_models(catalog, parameters)
reconstructionOptions = struct( ...
  "detection_threshold_deg", parameters.primary_grain_detection_deg, ...
  "min_grain_pixels", parameters.min_grain_pixels);
kernel = SO3DeLaValleePoussinKernel("halfwidth", ...
  parameters.texture_kernel_halfwidth_deg * degree);
models = cell(height(catalog), 2);
for scanIndex = 1:height(catalog)
  catalogRow = catalog(scanIndex,:);
  fprintf("ODF_RECOVERY sample=%s variant=%s\n", ...
    catalogRow.sample, catalogRow.variant);
  [ebsdFull, ~] = load_comprehensive_ebsd_scan(catalogRow);
  [grains, ebsdFull, ~] = reconstruct_comprehensive_grains( ...
    ebsdFull, reconstructionOptions);
  tiEbsd = ebsdFull("Ti-Hex");
  allTiGrains = grains("Ti-Hex");
  retainedGrain = allTiGrains.numPixel >= parameters.min_grain_pixels;
  tiGrains = allTiGrains(retainedGrain);
  assert(~isempty(tiEbsd) && ~isempty(tiGrains), ...
    "No retained Ti-Hex texture sources for %s %s.", ...
    catalogRow.sample, catalogRow.variant);
  grainAreas = abs(double(area(tiGrains)));
  grainAreas = grainAreas(:);
  assert(all(isfinite(grainAreas) & grainAreas > 0));

  pixelRbfOdf = calcDensity(tiEbsd.orientations, "kernel", kernel, ...
    "weights", ones(length(tiEbsd),1), "silent");
  grainRbfOdf = calcDensity(tiGrains.meanOrientation, ...
    "kernel", kernel, "weights", grainAreas / mean(grainAreas), ...
    "silent");
  models{scanIndex,1} = normalized_harmonic_odf( ...
    pixelRbfOdf, kernel.bandwidth);
  models{scanIndex,2} = normalized_harmonic_odf( ...
    grainRbfOdf, kernel.bandwidth);
  clear ebsdFull grains tiEbsd allTiGrains tiGrains
  clear pixelRbfOdf grainRbfOdf
end
end

function odf = normalized_harmonic_odf(inputOdf, targetBandwidth)
inputMean = double(mean(inputOdf));
assert(isscalar(inputMean) && isreal(inputMean) && ...
  isfinite(inputMean) && inputMean > 0);
normalizedInput = inputOdf / inputMean;
if isa(normalizedInput, "SO3FunHarmonic")
  odf = normalizedInput;
  if odf.bandwidth > targetBandwidth
    odf.bandwidth = targetBandwidth;
  end
else
  odf = SO3FunHarmonic(normalizedInput, "bandwidth", targetBandwidth);
end
assert(isa(odf, "SO3FunHarmonic") && ...
  odf.bandwidth <= targetBandwidth && ...
  abs(double(mean(odf)) - 1) < 1e-6);
end

function verify_rebuilt_texture_indices(models, catalog, textureDir, contract)
summaryPath = fullfile(textureDir, "texture_summary.csv");
assert(isfile(summaryPath), ...
  "ODF recovery requires existing texture summary: %s", summaryPath);
summary = readtable(summaryPath, "TextType", "string");
assert(isequal(string(summary.Properties.VariableNames), ...
  contract.summaryColumns.texture_summary));
assert(height(summary) == 2 * height(catalog));
weightingNames = ["pixel_weighted", "area_weighted_grain_mean"];
for scanIndex = 1:height(catalog)
  for weightingIndex = 1:2
    row = summary.sample == catalog.sample(scanIndex) & ...
      summary.variant == catalog.variant(scanIndex) & ...
      summary.weighting == weightingNames(weightingIndex);
    assert(nnz(row) == 1);
    rebuiltIndex = double(norm(models{scanIndex,weightingIndex})^2);
    registeredIndex = double(summary.texture_index(row));
    tolerance = 1e-8 * max(1, abs(registeredIndex));
    assert(abs(rebuiltIndex - registeredIndex) <= tolerance, ...
      ['Rebuilt ODF differs from registered texture index for ' ...
      '%s %s %s: rebuilt=%.15g registered=%.15g.'], ...
      char(catalog.sample(scanIndex)), char(catalog.variant(scanIndex)), ...
      char(weightingNames(weightingIndex)), rebuiltIndex, registeredIndex);
  end
end
end

function colorLimits = read_registered_texture_color_limits( ...
  textureDir, expectedScanCount, contract)
colorPath = fullfile(textureDir, "texture_plot_color_limits.csv");
assert(isfile(colorPath), ...
  "ODF recovery requires registered color limits: %s", colorPath);
colorLimits = readtable(colorPath, "TextType", "string");
assert(isequal(string(colorLimits.Properties.VariableNames), ...
  contract.summaryColumns.texture_plot_color_limits));
keys = colorLimits.plot_kind + "|" + colorLimits.weighting;
assert(height(colorLimits) == 6 && numel(unique(keys)) == 6);
assert(all(colorLimits.scan_count == expectedScanCount));
assert(all(colorLimits.grid_resolution_deg == ...
  contract.parameters.texture_grid_resolution_deg));
assert(all(isfinite(colorLimits.color_limit_min_mrd) & ...
  isfinite(colorLimits.color_limit_max_mrd) & ...
  colorLimits.color_limit_min_mrd == 0 & ...
  colorLimits.color_limit_max_mrd > 0));
for weighting = ["pixel_weighted", "area_weighted_grain_mean"]
  assert(nnz(colorLimits.plot_kind == "odf" & ...
    colorLimits.weighting == weighting) == 1);
end
end

function name = recovery_output_name(options)
name = "odf_sections.recovered.png";
if isfield(options, "texture_render_output_name")
  name = string(options.texture_render_output_name);
end
assert(isscalar(name) && strlength(name) > 0 && ...
  ~any(contains(name, ["/", "\"])), ...
  "texture_render_output_name must be a plain file name.");
[folder, base, extension] = fileparts(char(name));
assert(isempty(folder) && ~isempty(base) && strcmpi(extension, ".png"), ...
  "texture_render_output_name must be a plain PNG file name.");
end

function rows = selected_catalog_rows(options, catalogHeight)
if isfield(options, "catalog_rows")
  rows = double(options.catalog_rows(:));
else
  rows = (1:catalogHeight)';
end
assert(~isempty(rows) && all(isfinite(rows)) && ...
  all(rows == fix(rows)) && all(rows >= 1 & rows <= catalogHeight) && ...
  numel(unique(rows)) == numel(rows), ...
  "catalog_rows must be unique valid catalog row indices.");
end

function path = ensure_directory(root, name)
path = fullfile(root, name);
if ~isfolder(path)
  mkdir(path);
end
end

function meta = metadata_struct(catalogRow)
meta = struct("sample", string(catalogRow.sample), ...
  "diameter_mm", double(catalogRow.diameter_mm), ...
  "cold_reduction_percent", ...
  double(catalogRow.cold_reduction_percent), ...
  "variant", string(catalogRow.variant));
end

function output = add_scan_columns(input, catalogRow)
nRows = height(input);
output = addvars(input, repmat(catalogRow.sample, nRows, 1), ...
  repmat(catalogRow.diameter_mm, nRows, 1), ...
  repmat(catalogRow.cold_reduction_percent, nRows, 1), ...
  repmat(catalogRow.variant, nRows, 1), 'Before', 1, ...
  'NewVariableNames', {'sample', 'diameter_mm', ...
  'cold_reduction_percent', 'variant'});
end

function write_scan_table(scanTable, outputPath, firstScan)
if firstScan
  writetable(scanTable, outputPath);
else
  writetable(scanTable, outputPath, "WriteMode", "append");
end
end

function output = append_rows(output, rows)
if width(output) == 0
  output = rows;
else
  output = [output; rows];
end
end

function result = raw_denoised_distances(catalog, models, parameters)
sample = strings(0,1);
diameter_mm = zeros(0,1);
cold_reduction_percent = zeros(0,1);
weighting = strings(0,1);
kernel_halfwidth_deg = zeros(0,1);
harmonic_bandwidth = zeros(0,1);
odf_harmonic_l2_distance = zeros(0,1);
result = table(sample, diameter_mm, cold_reduction_percent, weighting, ...
  kernel_halfwidth_deg, harmonic_bandwidth, odf_harmonic_l2_distance);
weightingNames = ["pixel_weighted", "area_weighted_grain_mean"];
for sampleName = unique(catalog.sample, "stable")'
  rawRow = find(catalog.sample == sampleName & catalog.variant == "raw");
  denoisedRow = find(catalog.sample == sampleName & ...
    catalog.variant == "denoised");
  if numel(rawRow) ~= 1 || numel(denoisedRow) ~= 1
    continue
  end
  for weightingIndex = 1:2
    rawOdf = models{rawRow, weightingIndex};
    denoisedOdf = models{denoisedRow, weightingIndex};
    assert(isa(rawOdf, "SO3FunHarmonic") && ...
      isa(denoisedOdf, "SO3FunHarmonic"));
    assert(rawOdf.bandwidth == denoisedOdf.bandwidth);
    distance = double(norm(rawOdf - denoisedOdf));
    row = table(sampleName, double(catalog.diameter_mm(rawRow)), ...
      double(catalog.cold_reduction_percent(rawRow)), ...
      weightingNames(weightingIndex), ...
      double(parameters.texture_kernel_halfwidth_deg), ...
      double(rawOdf.bandwidth), distance, 'VariableNames', ...
      result.Properties.VariableNames);
    result = [result; row]; %#ok<AGROW>
  end
end
end

function plot_intragranular_trends(summary, outputFile)
selectors = { ...
  summary.kam_order == 1 & summary.kam_threshold_deg == 5, ...
  summary.kam_order == 2 & summary.kam_threshold_deg == 5, ...
  summary.kam_order == 1 & summary.kam_threshold_deg == 5, ...
  summary.kam_order == 1 & summary.kam_threshold_deg == 5};
metrics = ["kam_mean_deg", "kam_mean_deg", ...
  "grod_mean_deg", "gos_area_weighted_mean_deg"];
labels = ["KAM mean (deg), order 1 / 5 deg", ...
  "KAM mean (deg), order 2 / 5 deg", "GROD mean (deg)", ...
  "Area-weighted GOS mean (deg)"];
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [100 100 1400 950]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 2, 2, "Padding", "compact", ...
  "TileSpacing", "compact");
for metricIndex = 1:numel(metrics)
  axesHandle = nexttile(layout);
  hold(axesHandle, "on");
  for variant = ["raw", "denoised"]
    rows = selectors{metricIndex} & summary.variant == variant;
    if ~any(rows)
      continue
    end
    [xValues, order] = sort(summary.cold_reduction_percent(rows));
    yValues = summary.(metrics(metricIndex));
    yValues = yValues(rows);
    lineStyle = "-";
    if variant == "denoised"
      lineStyle = "--";
    end
    plot(axesHandle, xValues, yValues(order), "LineStyle", lineStyle, ...
      "Marker", "o", "LineWidth", 1.5, "DisplayName", variant);
  end
  xlabel(axesHandle, "Cold reduction (%)");
  ylabel(axesHandle, labels(metricIndex));
  grid(axesHandle, "on");
  legend(axesHandle, "Location", "best");
end
title(layout, ["Intragranular orientation-gradient/spread proxies", ...
  " (not direct dislocation-density measurements)"]);
exportgraphics(figureHandle, char(outputFile), "Resolution", 240, ...
  "BackgroundColor", "white");
clear cleanupFigure
end

function plot_texture_trends(summary, outputFile)
metrics = ["max_mrd", "texture_index", "m_index", ...
  "clipped_grid_entropy"];
labels = ["Maximum MRD", "Texture J-index", "M-index", ...
  "Clipped grid entropy"];
assert(all(ismember(metrics, string(summary.Properties.VariableNames))), ...
  "Texture trend table is missing a requested metric field.");
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [100 100 1400 950]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 2, 2, "Padding", "compact", ...
  "TileSpacing", "compact");
colors = lines(2);
weightingNames = ["pixel_weighted", "area_weighted_grain_mean"];
for metricIndex = 1:numel(metrics)
  axesHandle = nexttile(layout);
  hold(axesHandle, "on");
  for weightingIndex = 1:2
    for variant = ["raw", "denoised"]
      rows = summary.weighting == weightingNames(weightingIndex) & ...
        summary.variant == variant;
      if ~any(rows)
        continue
      end
      [xValues, order] = sort(summary.cold_reduction_percent(rows));
      yValues = summary.(metrics(metricIndex));
      yValues = yValues(rows);
      lineStyle = "-";
      if variant == "denoised"
        lineStyle = "--";
      end
      plot(axesHandle, xValues, yValues(order), ...
        "Color", colors(weightingIndex,:), "LineStyle", lineStyle, ...
        "Marker", "o", "LineWidth", 1.3, "DisplayName", ...
        weightingNames(weightingIndex) + " " + variant);
    end
  end
  xlabel(axesHandle, "Cold reduction (%)");
  ylabel(axesHandle, labels(metricIndex));
  grid(axesHandle, "on");
  legend(axesHandle, "Location", "best", "FontSize", 7, ...
    "Interpreter", "none");
end
title(layout, "Alpha-Ti texture metrics (raw primary; denoised sensitivity)");
exportgraphics(figureHandle, char(outputFile), "Resolution", 240, ...
  "BackgroundColor", "white");
clear cleanupFigure
end

function plot_texture_montage(models, catalog, colorLimits, plotKind, ...
  outputFile, parameters)
tempRoot = string(tempname);
mkdir(tempRoot);
cleanupTemp = onCleanup(@() remove_temp_directory(tempRoot));
weightingNames = ["pixel_weighted", "area_weighted_grain_mean"];
colorMaximums = zeros(1,2);
for weightingIndex = 1:2
  row = colorLimits.plot_kind == plotKind & ...
    colorLimits.weighting == weightingNames(weightingIndex);
  assert(nnz(row) == 1);
  colorMaximums(weightingIndex) = colorLimits.color_limit_max_mrd(row);
end
panelPaths = strings(height(catalog), 1);
for scanIndex = 1:height(catalog)
  panelPaths(scanIndex) = fullfile(tempRoot, ...
    sprintf("%02d_%s.png", scanIndex, plotKind));
  render_texture_panel(models(scanIndex,:), catalog(scanIndex,:), ...
    plotKind, panelPaths(scanIndex), colorMaximums, parameters);
end

columnCount = min(2, height(catalog));
rowCount = ceil(height(catalog) / columnCount);
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [20 20 2400 950 * rowCount]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, rowCount, columnCount, ...
  "Padding", "compact", "TileSpacing", "compact");
for scanIndex = 1:height(catalog)
  axesHandle = nexttile(layout);
  image(axesHandle, imread(panelPaths(scanIndex)));
  axis(axesHandle, "image");
  axis(axesHandle, "off");
end
exportgraphics(figureHandle, char(outputFile), "Resolution", 160, ...
  "BackgroundColor", "white");
clear cleanupFigure cleanupTemp
end

function render_texture_panel(models, catalogRow, plotKind, outputFile, ...
  colorMaximums, parameters)
if plotKind == "odf"
  render_odf_panel(models, catalogRow, outputFile, colorMaximums, ...
    parameters);
  return
end
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [50 50 1500 850]);
cleanupFigure = onCleanup(@() close(figureHandle));
  weightingLabels = ["Pixel weighted", "Area-wtd grain mean"];
odf = models{1};
cs = odf.CS;
switch plotKind
  case "pole"
    plotObjects = [Miller(0,0,0,1,cs), ...
      Miller(1,0,-1,0,cs), Miller(1,1,-2,0,cs)];
    columnLabels = ["{0001}", "{10-10}", "{11-20}"];
  case "ipf"
    plotObjects = [xvector, yvector, zvector];
    columnLabels = ["AD / X", "TD-RD / Y", "ND / Z"];
  otherwise
    error("Unknown texture plot kind: %s", plotKind);
end

for weightingIndex = 1:2
  for columnIndex = 1:3
    axesHandle = subplot(2, 3, ...
      (weightingIndex - 1) * 3 + columnIndex, "Parent", figureHandle);
    switch plotKind
      case "pole"
        plotPDF(models{weightingIndex}, plotObjects(columnIndex), ...
          "antipodal", "contourf", "silent", "parent", axesHandle, ...
          "resolution", parameters.texture_grid_resolution_deg * degree, ...
          "colorRange", [0 colorMaximums(weightingIndex)]);
      case "ipf"
        plotIPDF(models{weightingIndex}, plotObjects(columnIndex), ...
          "antipodal", "contourf", "silent", "parent", axesHandle, ...
          "resolution", parameters.texture_grid_resolution_deg * degree, ...
          "colorRange", [0 colorMaximums(weightingIndex)]);
    end
    title(axesHandle, weightingLabels(weightingIndex) + " | " + ...
      columnLabels(columnIndex), "Interpreter", "none", "FontSize", 9);
    if columnIndex == 3
      colorbarHandle = colorbar(axesHandle, "eastoutside");
      colorbarHandle.Label.String = sprintf("MRD (0-%.2f)", ...
        colorMaximums(weightingIndex));
    end
  end
end
sgtitle(figureHandle, sprintf( ...
  "%s %s | %.2f%% cold reduction | 5 deg kernel/grid | MRD", ...
  catalogRow.sample, catalogRow.variant, ...
  catalogRow.cold_reduction_percent), "Interpreter", "none");
drawnow;
exportgraphics(figureHandle, char(outputFile), "Resolution", 140, ...
  "BackgroundColor", "white");
clear cleanupFigure
end

function render_odf_panel(models, catalogRow, outputFile, ...
  colorMaximums, parameters)
tempRoot = string(tempname);
mkdir(tempRoot);
cleanupTemp = onCleanup(@() remove_temp_directory(tempRoot));
weightingLabels = ["Pixel weighted", "Area-weighted grain mean"];
rowPaths = strings(2,1);
for weightingIndex = 1:2
  rowPaths(weightingIndex) = fullfile(tempRoot, ...
    sprintf("odf_row_%d.png", weightingIndex));
  rowFigure = figure("Visible", "off", "Color", "white", ...
    "Position", [50 50 1500 520]);
  cleanupRow = onCleanup(@() close(rowFigure));
  plotSection(models{weightingIndex}, "phi2", [0 30 60] * degree, ...
    "contourf", "silent", "layout", [1 3], ...
    "resolution", parameters.texture_grid_resolution_deg * degree, ...
    "colorRange", [0 colorMaximums(weightingIndex)]);
  mtexColorbar("title", sprintf("MRD (0-%.2f)", ...
    colorMaximums(weightingIndex)));
  drawnow;
  exportgraphics(rowFigure, char(rowPaths(weightingIndex)), ...
    "Resolution", 140, "BackgroundColor", "white");
  clear cleanupRow
  crop_white_margins(rowPaths(weightingIndex), 12);
end

rowImages = cellfun(@imread, cellstr(rowPaths), ...
  "UniformOutput", false);
rowWidths = cellfun(@(x) size(x,2), rowImages);
rowHeights = cellfun(@(x) size(x,1), rowImages);
figureHeight = max(500, round(1500 * ...
  sum(rowHeights ./ rowWidths) + 100));
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [50 50 1500 figureHeight]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 2, 1, "Padding", "compact", ...
  "TileSpacing", "compact");
for weightingIndex = 1:2
  axesHandle = nexttile(layout);
  image(axesHandle, rowImages{weightingIndex});
  axis(axesHandle, "image");
  axis(axesHandle, "off");
  title(axesHandle, weightingLabels(weightingIndex), ...
    "Interpreter", "none", "FontSize", 10);
end
title(layout, sprintf( ...
  ['%s %s | %.2f%% cold reduction | 5 deg kernel/grid | ' ...
  'columns: phi2 = 0, 30, 60 deg | MRD'], ...
  char(catalogRow.sample), char(catalogRow.variant), ...
  catalogRow.cold_reduction_percent), "Interpreter", "none");
exportgraphics(figureHandle, char(outputFile), "Resolution", 140, ...
  "BackgroundColor", "white");
clear cleanupFigure cleanupTemp
end

function maxima = texture_color_maxima(models, plotKind, parameters)
maxima = zeros(1,2);
resolution = parameters.texture_grid_resolution_deg * degree;
for weightingIndex = 1:2
  for scanIndex = 1:size(models,1)
    odf = models{scanIndex, weightingIndex};
    switch plotKind
      case "pole"
        plotObjects = [Miller(0,0,0,1,odf.CS), ...
          Miller(1,0,-1,0,odf.CS), Miller(1,1,-2,0,odf.CS)];
        for objectIndex = 1:numel(plotObjects)
          density = calcPDF(odf, plotObjects(objectIndex), [], ...
            "antipodal");
          maxima(weightingIndex) = max(maxima(weightingIndex), ...
            double(max(density, "resolution", resolution)));
        end
      case "ipf"
        plotObjects = [xvector, yvector, zvector];
        for objectIndex = 1:numel(plotObjects)
          density = calcPDF(odf, [], plotObjects(objectIndex), ...
            "antipodal");
          maxima(weightingIndex) = max(maxima(weightingIndex), ...
            double(max(density, "resolution", resolution)));
        end
      case "odf"
        maxima(weightingIndex) = max(maxima(weightingIndex), ...
          double(max(odf, "resolution", resolution)));
      otherwise
        error("Unknown texture plot kind: %s", plotKind);
    end
  end
end
assert(all(isfinite(maxima) & maxima > 0), ...
  "Texture color scale contains an invalid maximum.");
end

function result = texture_plot_color_limits(models, parameters)
plot_kind = strings(0,1);
weighting = strings(0,1);
scan_count = zeros(0,1);
grid_resolution_deg = zeros(0,1);
color_limit_min_mrd = zeros(0,1);
color_limit_max_mrd = zeros(0,1);
result = table(plot_kind, weighting, scan_count, grid_resolution_deg, ...
  color_limit_min_mrd, color_limit_max_mrd);
weightingNames = ["pixel_weighted", "area_weighted_grain_mean"];
for plotKind = ["pole", "ipf", "odf"]
  maxima = texture_color_maxima(models, plotKind, parameters);
  rows = table(repmat(plotKind,2,1), weightingNames', ...
    repmat(size(models,1),2,1), ...
    repmat(parameters.texture_grid_resolution_deg,2,1), ...
    zeros(2,1), maxima', ...
    'VariableNames', result.Properties.VariableNames);
  result = [result; rows]; %#ok<AGROW>
end
keys = result.plot_kind + "|" + result.weighting;
assert(height(result) == 6 && numel(unique(keys)) == 6);
end

function crop_white_margins(path, padding)
imageData = imread(path);
if ndims(imageData) == 2
  nonwhite = imageData < 250;
else
  nonwhite = any(imageData < 250, 3);
end
[rowIndex, columnIndex] = find(nonwhite);
assert(~isempty(rowIndex), "Texture panel is blank: %s", path);
firstRow = max(1, min(rowIndex) - padding);
lastRow = min(size(imageData,1), max(rowIndex) + padding);
firstColumn = max(1, min(columnIndex) - padding);
lastColumn = min(size(imageData,2), max(columnIndex) + padding);
imwrite(imageData(firstRow:lastRow, firstColumn:lastColumn, :), path);
end

function remove_temp_directory(path)
if isfolder(path)
  rmdir(path, "s");
end
end
