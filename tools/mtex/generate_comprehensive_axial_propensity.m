function [byGrain, summary] = generate_comprehensive_axial_propensity( ...
  scanRoot, outputRoot, options)
%GENERATE_COMPREHENSIVE_AXIAL_PROPENSITY Export AD tensile propensity data.

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
artifactDir = fullfile(outputRoot, "06_axial_propensity");
if ~isfolder(artifactDir)
  mkdir(artifactDir);
end

reconstructionOptions = struct( ...
  "detection_threshold_deg", ...
  contract.parameters.primary_grain_detection_deg, ...
  "min_grain_pixels", contract.parameters.min_grain_pixels);
byGrain = table();
summary = table();

for scanIndex = 1:height(catalog)
  catalogRow = catalog(scanIndex, :);
  fprintf("AXIAL_PROPENSITY sample=%s variant=%s\n", ...
    catalogRow.sample, catalogRow.variant);
  [ebsdFull, ~] = load_comprehensive_ebsd_scan(catalogRow);
  [grains, ~, ~] = reconstruct_comprehensive_grains( ...
    ebsdFull, reconstructionOptions);
  retained = grains("Ti-Hex");
  retained = retained(retained.numPixel >= ...
    reconstructionOptions.min_grain_pixels);
  assert(~isempty(retained), ...
    "No retained Ti-Hex grains for %s %s.", ...
    catalogRow.sample, catalogRow.variant);

  [scanByGrain, scanSummary] = compute_axial_propensity_metrics( ...
    retained.meanOrientation, double(retained.id(:)), ...
    double(area(retained)), options);
  scanByGrain = add_scan_columns(scanByGrain, catalogRow);
  scanSummary = add_scan_columns(scanSummary, catalogRow);
  assert(isequal(string(scanSummary.Properties.VariableNames), ...
    contract.summaryColumns.axial_propensity_summary));
  if scanIndex == 1
    byGrain = scanByGrain;
    summary = scanSummary;
  else
    byGrain = [byGrain; scanByGrain]; %#ok<AGROW>
    summary = [summary; scanSummary]; %#ok<AGROW>
  end
  clear ebsdFull grains retained scanByGrain scanSummary
end

assert(height(summary) > 0 && height(byGrain) > 0);
assert(all(byGrain.max_abs_schmid >= 0 & ...
  byGrain.max_abs_schmid <= 0.5 + 1e-10));
writetable(byGrain, fullfile(artifactDir, ...
  "axial_propensity_by_grain.csv"));
writetable(summary, fullfile(artifactDir, ...
  "axial_propensity_summary.csv"));
plot_propensity_trends(summary, fullfile(artifactDir, ...
  "axial_propensity_trends.png"));
end

function output = add_scan_columns(input, catalogRow)
nRows = height(input);
output = addvars(input, ...
  repmat(catalogRow.sample, nRows, 1), ...
  repmat(catalogRow.diameter_mm, nRows, 1), ...
  repmat(catalogRow.cold_reduction_percent, nRows, 1), ...
  repmat(catalogRow.variant, nRows, 1), ...
  'Before', 1, 'NewVariableNames', ...
  {'sample', 'diameter_mm', 'cold_reduction_percent', 'variant'});
end

function plot_propensity_trends(summary, outputFile)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [100 100 1400 900]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 2, 1, "Padding", "compact", ...
  "TileSpacing", "compact");
families = unique(summary.family, "stable");
colors = lines(numel(families));
lineStyles = struct("raw", "-", "denoised", "--");

metricNames = ["area_weighted_mean_max_schmid"; ...
  "area_fraction_schmid_ge_0_4"];
yLabels = ["Area-weighted mean max |m|"; ...
  "Area fraction with max |m| >= 0.4"];
for metricIndex = 1:numel(metricNames)
  axesHandle = nexttile(layout);
  hold(axesHandle, "on");
  for familyIndex = 1:numel(families)
    for variant = ["raw", "denoised"]
      rows = summary.family == families(familyIndex) & ...
        summary.variant == variant;
      [x, order] = sort(summary.cold_reduction_percent(rows));
      values = summary.(metricNames(metricIndex));
      values = values(rows);
      displayName = families(familyIndex) + " " + variant;
      plot(axesHandle, x, values(order), ...
        "LineStyle", lineStyles.(variant), ...
        "Color", colors(familyIndex, :), "Marker", "o", ...
        "LineWidth", 1.4, "MarkerSize", 4, ...
        "DisplayName", displayName);
    end
  end
  xlabel(axesHandle, "Cold reduction (%)");
  ylabel(axesHandle, yLabels(metricIndex));
  grid(axesHandle, "on");
  if metricIndex == 1
    title(axesHandle, ...
      "AD uniaxial-tension geometric propensity (raw primary; denoised sensitivity)");
    legend(axesHandle, "Location", "eastoutside", ...
      "Interpreter", "none");
  end
end

drawnow;
exportgraphics(figureHandle, char(outputFile), "Resolution", 300, ...
  "BackgroundColor", "white");
renderedImage = imread(outputFile);
imwrite(renderedImage, outputFile, "png");
clear cleanupFigure
end
