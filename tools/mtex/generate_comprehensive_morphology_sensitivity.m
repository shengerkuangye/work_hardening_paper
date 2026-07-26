function summary = generate_comprehensive_morphology_sensitivity( ...
  scanRoot, outputRoot)
%GENERATE_COMPREHENSIVE_MORPHOLOGY_SENSITIVITY Raw-only shape sensitivity.

arguments
  scanRoot (1,1) string
  outputRoot (1,1) string
end

assert(isfolder(scanRoot), "EBSD scan folder not found: %s", scanRoot);
assert(~isempty(which("EBSD")), ...
  "MTEX is not loaded in the current MATLAB session.");
formalFile = fullfile(outputRoot, "02_grain_morphology", ...
  "grain_morphology_summary.csv");
assert(isfile(formalFile), ...
  "Formal morphology summary not found: %s", formalFile);

catalog = comprehensive_ebsd_catalog(scanRoot);
catalog = catalog(catalog.variant == "raw", :);
assert(height(catalog) == 6);
sensitivityDir = fullfile(outputRoot, "09_sensitivity");
if ~isfolder(sensitivityDir)
  mkdir(sensitivityDir);
end

domainDetectionDeg = 2;
domainMinPixels = [1 3 5 10];
hagbDetectionDeg = 15;
hagbMinPixels = 5;
summary = table();

for scanIndex = 1:height(catalog)
  catalogRow = catalog(scanIndex,:);
  fprintf("MORPHOLOGY_SENSITIVITY sample=%s variant=raw\n", ...
    catalogRow.sample);
  [ebsdFull, ~] = load_comprehensive_ebsd_scan(catalogRow);
  meta = catalogRow(:, {'sample','diameter_mm', ...
    'cold_reduction_percent','variant'});

  domainOptions = struct("detection_threshold_deg", ...
    domainDetectionDeg, "min_grain_pixels", 1);
  [domainObjects, ~] = reconstruct_comprehensive_grains( ...
    ebsdFull, domainOptions);
  domainRows = compute_morphology_sensitivity_rows( ...
    domainObjects("Ti-Hex"), meta, "orientation_domain", ...
    domainDetectionDeg, domainMinPixels);
  summary = append_rows(summary, domainRows);
  clear domainObjects domainRows

  hagbOptions = struct("detection_threshold_deg", ...
    hagbDetectionDeg, "min_grain_pixels", hagbMinPixels);
  [hagbObjects, ~] = reconstruct_comprehensive_grains( ...
    ebsdFull, hagbOptions);
  hagbRows = compute_morphology_sensitivity_rows( ...
    hagbObjects("Ti-Hex"), meta, "hagb_grain", ...
    hagbDetectionDeg, hagbMinPixels);
  summary = append_rows(summary, hagbRows);
  clear ebsdFull hagbObjects hagbRows
end

assert(height(summary) == 30);
assert(all(summary.variant == "raw"));
validate_against_formal(summary, formalFile);

csvFile = fullfile(sensitivityDir, ...
  "morphology_sensitivity_summary.csv");
figureFile = fullfile(sensitivityDir, ...
  "morphology_sensitivity_trends.png");
writetable(summary, csvFile);
plot_sensitivity(summary, figureFile);
fprintf("MORPHOLOGY_SENSITIVITY_SUMMARY=%s\n", csvFile);
fprintf("MORPHOLOGY_SENSITIVITY_FIGURE=%s\n", figureFile);
end

function output = append_rows(output, rows)
if width(output) == 0
  output = rows;
else
  output = [output; rows];
end
end

function validate_against_formal(summary, formalFile)
formal = readtable(formalFile, 'TextType', 'string');
formal = formal(formal.variant == "raw", :);
primary = summary(summary.object_type == "orientation_domain" & ...
  summary.grain_detection_deg == 2 & ...
  summary.min_grain_pixels == 5, :);
primary = sortrows(primary, "cold_reduction_percent");
formal = sortrows(formal, "cold_reduction_percent");
assert(height(primary) == 6 && height(formal) == 6);
assert(isequal(primary.sample, formal.sample));
assert(isequal(primary.retained_object_count, formal.grain_count));

differences = [ ...
  primary.retained_area_um2 - formal.total_area_um2, ...
  primary.ecd_number_median_um - formal.ecd_number_median_um, ...
  primary.ecd_area_weighted_median_um - ...
    formal.ecd_area_weighted_median_um, ...
  primary.aspect_ratio_number_median - ...
    formal.aspect_ratio_number_median, ...
  primary.aspect_ratio_area_weighted_median - ...
    formal.aspect_ratio_area_weighted_median, ...
  primary.long_axis_ad_angle_area_weighted_median_deg - ...
    formal.long_axis_ad_angle_area_weighted_median_deg];
assert(all(isfinite(differences), "all"));
assert(max(abs(differences), [], "all") < 1e-10, ...
  "2 degree/minPixel 5 rows differ from formal morphology summary.");
end

function plot_sensitivity(summary, outputFile)
domainRows = summary(summary.object_type == "orientation_domain", :);
hagbRows = summary(summary.object_type == "hagb_grain", :);
minPixels = unique(domainRows.min_grain_pixels, "stable");
colors = lines(numel(minPixels));
metricNames = ["retained_area_fraction", ...
  "ecd_area_weighted_median_um", ...
  "aspect_ratio_area_weighted_median", ...
  "long_axis_ad_angle_area_weighted_median_deg", ...
  "retained_object_count"];
yLabels = ["Retained area fraction", ...
  "Area-weighted median ECD (um)", ...
  "Area-weighted median aspect ratio", ...
  "Area-weighted median long-axis--AD angle (deg)", ...
  "Retained domain / grain count"];

figureHandle = figure('Visible', 'off', ...
  'Position', [100 100 1700 950], 'Color', 'w');
layout = tiledlayout(figureHandle, 2, 3, ...
  'TileSpacing', 'compact', 'Padding', 'compact');
legendHandles = gobjects(numel(minPixels) + 1, 1);
legendLabels = strings(numel(minPixels) + 1, 1);

for metricIndex = 1:numel(metricNames)
  axesHandle = nexttile(layout);
  hold(axesHandle, 'on');
  for minIndex = 1:numel(minPixels)
    rows = domainRows(domainRows.min_grain_pixels == ...
      minPixels(minIndex), :);
    rows = sortrows(rows, "cold_reduction_percent");
    lineHandle = plot(axesHandle, rows.cold_reduction_percent, ...
      rows.(metricNames(metricIndex)), '-o', 'LineWidth', 1.6, ...
      'MarkerSize', 5, 'Color', colors(minIndex,:));
    if metricIndex == 1
      legendHandles(minIndex) = lineHandle;
      legendLabels(minIndex) = sprintf( ...
        "2 deg domains, minPixel=%d", minPixels(minIndex));
    end
  end
  hagbRows = sortrows(hagbRows, "cold_reduction_percent");
  lineHandle = plot(axesHandle, hagbRows.cold_reduction_percent, ...
    hagbRows.(metricNames(metricIndex)), '--s', 'LineWidth', 1.8, ...
    'MarkerSize', 5, 'Color', [0 0 0]);
  if metricIndex == 1
    legendHandles(end) = lineHandle;
    legendLabels(end) = "15 deg HAGB grains, minPixel=5";
  end
  xlabel(axesHandle, "Cold reduction (%)");
  ylabel(axesHandle, yLabels(metricIndex));
  grid(axesHandle, 'on');
  box(axesHandle, 'on');
  if metricNames(metricIndex) == "retained_object_count"
    set(axesHandle, 'YScale', 'log');
  end
end

legendAxes = nexttile(layout);
axis(legendAxes, 'off');
legend(legendAxes, legendHandles, legendLabels, ...
  'Location', 'northwest');
title(layout, ["Raw EBSD morphology sensitivity: " ...
  "2 deg orientation domains and 15 deg HAGB grains"]);
exportgraphics(figureHandle, outputFile, 'Resolution', 300);
close(figureHandle);
end
