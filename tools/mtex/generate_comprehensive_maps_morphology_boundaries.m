function [morphology, morphologySummary, boundarySegments, ...
  boundarySummary, boundaryDistribution, sensitivity] = ...
  generate_comprehensive_maps_morphology_boundaries(scanRoot, ...
  outputRoot)
%GENERATE_COMPREHENSIVE_MAPS_MORPHOLOGY_BOUNDARIES Export Task 3 results.

arguments
  scanRoot (1,1) string
  outputRoot (1,1) string
end

assert(isfolder(scanRoot), "EBSD scan folder not found: %s", scanRoot);
assert(~isempty(which("EBSD")), ...
  "MTEX is not loaded in the current MATLAB session.");
contract = comprehensive_ebsd_output_contract();
catalog = comprehensive_ebsd_catalog(scanRoot);
mapDir = ensure_directory(outputRoot, "01_standard_maps");
morphologyDir = ensure_directory(outputRoot, "02_grain_morphology");
boundaryDir = ensure_directory(outputRoot, "03_boundaries");

parameters = contract.parameters;
primaryFloorDeg = parameters.primary_grain_detection_deg;
classificationDeg = parameters.boundary_classification_deg;
minGrainPixels = parameters.min_grain_pixels;
detectionFloorsDeg = parameters.detection_sensitivity_deg;
reconstructionOptions = struct( ...
  "detection_threshold_deg", primaryFloorDeg, ...
  "min_grain_pixels", minGrainPixels);

morphology = table();
morphologySummary = table();
morphologyQuantiles = table();
boundarySegments = table();
boundarySummary = table();
boundaryDistribution = table();
sensitivity = table();

for scanIndex = 1:height(catalog)
  catalogRow = catalog(scanIndex, :);
  fprintf("MAPS_MORPHOLOGY_BOUNDARIES sample=%s variant=%s\n", ...
    catalogRow.sample, catalogRow.variant);
  [ebsdFull, scanMeta] = load_comprehensive_ebsd_scan(catalogRow);
  [morphologyGrains, ebsdMorphology] = ...
    reconstruct_comprehensive_grains(ebsdFull, reconstructionOptions);
  tiGrains = morphologyGrains("Ti-Hex");
  metricMeta = common_metadata(catalogRow);
  morphologyOptions = struct( ...
    "grain_detection_deg", primaryFloorDeg, ...
    "min_grain_pixels", minGrainPixels);
  [scanMorphology, scanMorphologySummary, scanQuantiles] = ...
    compute_grain_morphology_metrics(tiGrains, metricMeta, ...
    morphologyOptions);
  scanQuantiles = add_common_columns(scanQuantiles, catalogRow);
  morphology = append_rows(morphology, scanMorphology);
  morphologySummary = append_rows(morphologySummary, ...
    scanMorphologySummary);
  morphologyQuantiles = append_rows(morphologyQuantiles, scanQuantiles);

  sourceFloorDeg = min(detectionFloorsDeg);
  [boundarySource, boundaryObjects] = boundary_source( ...
    ebsdFull, sourceFloorDeg, classificationDeg, minGrainPixels, ...
    scanMeta, catalogRow);
  indexedAreaUm2 = length(ebsdFull("Ti-Hex")) * ...
    scanMeta.x_step_um * scanMeta.y_step_um;
  for floorDeg = detectionFloorsDeg
    fprintf("  boundary detection floor=%.1f deg\n", floorDeg);
    boundaryOptions = struct( ...
      "grain_detection_deg", floorDeg, ...
      "detection_floor_deg", floorDeg, ...
      "classification_deg", classificationDeg, ...
      "indexed_area_um2", indexedAreaUm2, ...
      "min_boundary_axis_deg", parameters.min_boundary_axis_deg, ...
      "twin_candidate_tolerance_deg", ...
      parameters.twin_candidate_tolerance_deg);
    [scanSegments, scanBoundarySummary, scanDistribution] = ...
      compute_boundary_network_metrics(boundarySource, metricMeta, ...
      boundaryOptions);
    sensitivity = append_rows(sensitivity, scanBoundarySummary);
    if floorDeg == primaryFloorDeg
      boundarySegments = append_rows(boundarySegments, scanSegments);
      boundarySummary = append_rows(boundarySummary, ...
        scanBoundarySummary);
      boundaryDistribution = append_rows(boundaryDistribution, ...
        scanDistribution);
    end
    clear scanSegments scanBoundarySummary scanDistribution
  end

  mapPath = fullfile(mapDir, catalogRow.sample + "_" + ...
    catalogRow.variant + "_maps.png");
  export_standard_maps(ebsdMorphology, boundaryObjects, ...
    boundarySource, catalogRow, parameters, mapPath);
  clear ebsdFull ebsdMorphology morphologyGrains tiGrains
  clear boundarySource boundaryObjects
end

assert(height(morphologySummary) == height(catalog));
assert(height(boundarySummary) == height(catalog));
assert(height(sensitivity) == ...
  height(catalog) * numel(detectionFloorsDeg));
writetable(morphology, fullfile(morphologyDir, ...
  "grain_morphology_by_grain.csv"));
writetable(morphologySummary, fullfile(morphologyDir, ...
  "grain_morphology_summary.csv"));
writetable(morphologyQuantiles, fullfile(morphologyDir, ...
  "grain_morphology_quantiles.csv"));
writetable(boundarySegments, fullfile(boundaryDir, ...
  "boundary_segments.csv"));
writetable(boundarySummary, fullfile(boundaryDir, ...
  "boundary_summary.csv"));
writetable(boundaryDistribution, fullfile(boundaryDir, ...
  "boundary_angle_distributions.csv"));
writetable(sensitivity, fullfile(boundaryDir, ...
  "boundary_detection_sensitivity.csv"));
plot_morphology_trends(morphologySummary, fullfile(morphologyDir, ...
  "grain_morphology_trends.png"));
plot_boundary_trends(boundarySummary, fullfile(boundaryDir, ...
  "boundary_trends.png"));
end

function path = ensure_directory(root, name)
path = fullfile(root, name);
if ~isfolder(path)
  mkdir(path);
end
end

function meta = common_metadata(catalogRow)
meta = catalogRow(:, {'sample','diameter_mm', ...
  'cold_reduction_percent','variant'});
end

function output = add_common_columns(input, catalogRow)
nRows = height(input);
output = addvars(input, repmat(catalogRow.sample, nRows, 1), ...
  repmat(catalogRow.diameter_mm, nRows, 1), ...
  repmat(catalogRow.cold_reduction_percent, nRows, 1), ...
  repmat(catalogRow.variant, nRows, 1), 'Before', 1, ...
  'NewVariableNames', {'sample','diameter_mm', ...
  'cold_reduction_percent','variant'});
end

function output = append_rows(output, rows)
if width(output) == 0
  output = rows;
else
  output = [output; rows];
end
end

function [source, objects] = boundary_source(ebsdFull, sourceFloorDeg, ...
  classificationDeg, minGrainPixels, scanMeta, catalogRow)
[grains, grainId] = calcGrains(ebsdFull, 'unitCell', 'threshold', ...
  [classificationDeg sourceFloorDeg] * degree, ...
  'minPixel', minGrainPixels);
ebsdFull.grainId = grainId;
objects.inner = grains.innerBoundary("Ti-Hex", "Ti-Hex");
objects.outer = grains.boundary("Ti-Hex", "Ti-Hex");
assert(~isempty(objects.inner) && ~isempty(objects.outer), ...
  "No Ti-Hex boundary network for %s %s at %.1f degrees.", ...
  catalogRow.sample, catalogRow.variant, sourceFloorDeg);

source.inner_theta_deg = ...
  double(angle(objects.inner.misorientation) / degree);
source.inner_length_um = double(objects.inner.segLength(:));
source.inner_endpoint_ids = double(objects.inner.ebsdId);
source.outer_theta_deg = ...
  double(angle(objects.outer.misorientation) / degree);
source.outer_length_um = double(objects.outer.segLength(:));
source.outer_endpoint_ids = double(objects.outer.ebsdId);

source.inner_axis_xyz = boundary_axes(objects.inner);
source.outer_axis_xyz = boundary_axes(objects.outer);
tiCrystalSymmetry = ebsdFull("Ti-Hex").CS;
extensionTwin = orientation.map( ...
  Miller(1, -1, 0, 1, tiCrystalSymmetry), ...
  Miller(1, 0, -1, -1, tiCrystalSymmetry), ...
  Miller(0, 1, -1, 1, tiCrystalSymmetry, 'uvw'), ...
  Miller(1, -1, 0, 1, tiCrystalSymmetry, 'uvw'));
source.inner_twin_deviation_deg = double( ...
  angle(objects.inner.misorientation, extensionTwin) / degree);
source.outer_twin_deviation_deg = double( ...
  angle(objects.outer.misorientation, extensionTwin) / degree);

endpointIds = [source.inner_endpoint_ids; source.outer_endpoint_ids];
topologyAudit = audit_native_grid_pairs(endpointIds, ...
  double(ebsdFull.id(:)), ...
  [double(ebsdFull.x(:)), double(ebsdFull.y(:))], ...
  scanMeta.x_step_um, scanMeta.y_step_um, ...
  catalogRow.sample + " " + catalogRow.variant + ...
  " Ti-Hex/Ti-Hex source faces");
assert(topologyAudit.endpoint_pair_count == ...
  numel(source.inner_theta_deg) + numel(source.outer_theta_deg));
assert(topologyAudit.nonlocal_endpoint_pair_count == 0);
end

function values = boundary_axes(boundary)
directions = axis(boundary.misorientation);
values = double([directions.x(:), directions.y(:), directions.z(:)]);
end

function export_standard_maps(ebsdFull, boundaryObjects, source, ...
  catalogRow, parameters, outputFile)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [20 20 1900 1650]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 3, 3, "Padding", "compact", ...
  "TileSpacing", "compact");
tiEbsd = ebsdFull("Ti-Hex");
directions = {xvector, yvector, zvector};
directionNames = ["IPF-X / AD", "IPF-Y / TD-RD", "IPF-Z / ND"];
for directionIndex = 1:3
  colorKey = ipfHSVKey(tiEbsd);
  colorKey.inversePoleFigureDirection = directions{directionIndex};
  colors = colorKey.orientation2color(tiEbsd.orientations);
  axesHandle = nexttile(layout);
  plot_rgb_map(axesHandle, ebsdFull, tiEbsd, colors);
  title(axesHandle, directionNames(directionIndex));
end

qualityMaps = {double(ebsdFull.prop.bc(:)), ...
  double(ebsdFull.prop.mad(:)), double(ebsdFull.phaseId(:))};
qualityTitles = ["Band contrast", "MAD", "Phase / indexing"];
for mapIndex = 1:3
  axesHandle = nexttile(layout);
  plot_scalar_map(axesHandle, ebsdFull, qualityMaps{mapIndex});
  title(axesHandle, qualityTitles(mapIndex));
  colorbar(axesHandle);
end

axesHandle = nexttile(layout);
grainColors = categorical_grain_colors(double(ebsdFull.grainId(:)));
plot_rgb_map(axesHandle, ebsdFull, ebsdFull, grainColors);
title(axesHandle, "Grain ID (2 degree reconstruction; categorical)");

axesHandle = nexttile(layout);
plot_scalar_map(axesHandle, ebsdFull, double(ebsdFull.prop.bc(:)));
colormap(axesHandle, gray(256));
hold(axesHandle, "on");
innerTheta = source.inner_theta_deg;
outerTheta = source.outer_theta_deg;
draw_boundary_segments(axesHandle, ...
  boundaryObjects.inner(innerTheta >= 2 & innerTheta < 5), ...
  [0.00 0.75 1.00], 0.8);
draw_boundary_segments(axesHandle, ...
  boundaryObjects.inner(innerTheta >= 5 & innerTheta < 15), ...
  [1.00 0.65 0.00], 0.9);
draw_boundary_segments(axesHandle, ...
  boundaryObjects.outer(outerTheta >= 15), [0.85 0.05 0.05], 1.0);
twinMask = outerTheta >= parameters.min_boundary_axis_deg & ...
  source.outer_twin_deviation_deg <= ...
  parameters.twin_candidate_tolerance_deg;
draw_boundary_segments(axesHandle, boundaryObjects.outer(twinMask), ...
  [0.80 0.00 0.80], 1.8);
legendHandles = [ ...
  plot(axesHandle, nan, nan, "Color", [0.00 0.75 1.00], ...
    "LineWidth", 1.5), ...
  plot(axesHandle, nan, nan, "Color", [1.00 0.65 0.00], ...
    "LineWidth", 1.5), ...
  plot(axesHandle, nan, nan, "Color", [0.85 0.05 0.05], ...
    "LineWidth", 1.5), ...
  plot(axesHandle, nan, nan, "Color", [0.80 0.00 0.80], ...
    "LineWidth", 1.8)];
legend(axesHandle, legendHandles, ...
  ["2-<5 deg", "5-<15 deg", ">=15 deg", ...
  "extension-twin angular-axis candidate (non-unique)"], ...
  "Location", "southoutside", "FontSize", 7);
title(axesHandle, "BC with unsmoothed boundary classes");

axesHandle = nexttile(layout);
indexedMap = double(ebsdFull.isIndexed(:));
plot_scalar_map(axesHandle, ebsdFull, indexedMap);
colormap(axesHandle, [0.85 0.85 0.85; 0.15 0.45 0.85]);
clim(axesHandle, [0 1]);
title(axesHandle, "Indexing mask (0 unindexed, 1 indexed)");
colorbar(axesHandle);
title(layout, sprintf("%s %s | %.2f%% cold reduction | AD horizontal", ...
  catalogRow.sample, catalogRow.variant, ...
  catalogRow.cold_reduction_percent), "Interpreter", "none");
drawnow;
exportgraphics(figureHandle, char(outputFile), "Resolution", 180, ...
  "BackgroundColor", "white");
clear cleanupFigure
end

function plot_rgb_map(axesHandle, ebsdFull, subset, rgbValues)
[xValues, yValues, linearIndices] = native_grid_indices(ebsdFull);
rgbImage = ones(numel(yValues), numel(xValues), 3);
[found, fullRows] = ismember(double(subset.id(:)), ...
  double(ebsdFull.id(:)));
assert(all(found));
for channelIndex = 1:3
  channel = rgbImage(:, :, channelIndex);
  channel(linearIndices(fullRows)) = rgbValues(:, channelIndex);
  rgbImage(:, :, channelIndex) = channel;
end
imagesc(axesHandle, xValues, yValues, rgbImage);
format_map_axes(axesHandle);
end

function plot_scalar_map(axesHandle, ebsdFull, values)
[xValues, yValues, linearIndices] = native_grid_indices(ebsdFull);
assert(numel(values) == length(ebsdFull));
imageData = nan(numel(yValues), numel(xValues));
imageData(linearIndices) = values;
imagesc(axesHandle, xValues, yValues, imageData);
format_map_axes(axesHandle);
end

function colors = categorical_grain_colors(grainIds)
grainIds = double(grainIds(:));
valid = isfinite(grainIds) & grainIds > 0;
colors = repmat([0.78 0.78 0.78], numel(grainIds), 1);
goldenRatioConjugate = (sqrt(5) - 1) / 2;
hue = mod(grainIds(valid) * goldenRatioConjugate, 1);
saturation = 0.55 + 0.35 * mod(grainIds(valid) * 0.381966, 1);
value = 0.75 + 0.25 * mod(grainIds(valid) * 0.754878, 1);
colors(valid, :) = hsv2rgb([hue, saturation, value]);
end

function [xValues, yValues, linearIndices] = native_grid_indices(ebsd)
xCoordinates = double(ebsd.x(:));
yCoordinates = double(ebsd.y(:));
xValues = unique(xCoordinates);
yValues = unique(yCoordinates);
[xFound, xIndex] = ismember(xCoordinates, xValues);
[yFound, yIndex] = ismember(yCoordinates, yValues);
assert(all(xFound & yFound));
linearIndices = sub2ind([numel(yValues), numel(xValues)], ...
  yIndex, xIndex);
assert(numel(unique(linearIndices)) == length(ebsd));
end

function format_map_axes(axesHandle)
axis(axesHandle, "image");
set(axesHandle, "YDir", "normal", "FontSize", 7);
xlabel(axesHandle, "AD (um)");
ylabel(axesHandle, "TD/RD (um)");
end

function draw_boundary_segments(axesHandle, boundary, color, lineWidth)
if isempty(boundary)
  return
end
vertices = boundary.allV.xyz;
faces = boundary.F;
xCoordinates = reshape(vertices(faces.', 1), 2, []);
yCoordinates = reshape(vertices(faces.', 2), 2, []);
separator = nan(1, size(xCoordinates, 2));
xCoordinates = [xCoordinates; separator];
yCoordinates = [yCoordinates; separator];
line(axesHandle, xCoordinates(:), yCoordinates(:), ...
  "Color", color, "LineWidth", lineWidth);
end

function plot_morphology_trends(summary, outputFile)
metrics = ["ecd_area_weighted_median_um", ...
  "aspect_ratio_area_weighted_median", ...
  "long_axis_ad_angle_area_weighted_median_deg", ...
  "boundary_touching_grain_fraction"];
labels = ["Area-weighted median ECD (um)", ...
  "Area-weighted median aspect ratio", ...
  "Area-weighted median long-axis angle from AD (deg)", ...
  "Boundary-touching grain fraction"];
plot_trend_panels(summary, metrics, labels, ...
  "Unsmoothed Ti-Hex grain morphology", outputFile);
end

function plot_boundary_trends(summary, outputFile)
metrics = ["boundary_line_density_per_um", ...
  "lagb_2_15_length_fraction", "hagb_ge15_length_fraction", ...
  "twin_candidate_length_fraction"];
labels = ["Boundary line density (um/um^2)", ...
  "2-<15 deg length fraction", ">=15 deg length fraction", ...
  "Extension-twin candidate length fraction"];
plot_trend_panels(summary, metrics, labels, ...
  "Length-weighted Ti-Hex boundary network", outputFile);
end

function plot_trend_panels(summary, metrics, labels, figureTitle, ...
  outputFile)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [100 100 1400 950]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 2, 2, "Padding", "compact", ...
  "TileSpacing", "compact");
for metricIndex = 1:numel(metrics)
  axesHandle = nexttile(layout);
  hold(axesHandle, "on");
  for variant = ["raw", "denoised"]
    rows = summary.variant == variant;
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
title(layout, figureTitle + " (raw primary; denoised sensitivity)");
drawnow;
exportgraphics(figureHandle, char(outputFile), "Resolution", 240, ...
  "BackgroundColor", "white");
clear cleanupFigure
end
