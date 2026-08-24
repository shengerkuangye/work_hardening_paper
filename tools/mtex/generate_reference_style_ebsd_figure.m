function summary = generate_reference_style_ebsd_figure(scanRoot, outputDir)
%GENERATE_REFERENCE_STYLE_EBSD_FIGURE Recreate the visual logic of Fig. 6.
% Columns compare the undeformed and final cold-deformed states. Rows show
% IPF-AD, classified boundaries, and KAM maps. The right column contains
% paired grain-size, boundary-misorientation, and KAM distributions.

arguments
  scanRoot (1,1) string
  outputDir (1,1) string
end

assert(isfolder(scanRoot), "EBSD scan folder not found: %s", scanRoot);
assert(~isempty(which("EBSD")), "MTEX must be loaded before this function.");
if ~isfolder(outputDir)
  mkdir(outputDir);
end

parameters = struct( ...
  "boundary_detection_deg", 2, ...
  "hagb_threshold_deg", 15, ...
  "min_grain_pixels", 5, ...
  "kam_order", 1, ...
  "kam_cutoff_deg", 5, ...
  "grain_bin_width_um", 2, ...
  "misorientation_bin_width_deg", 2, ...
  "kam_bin_width_deg", 0.2);

catalog = comprehensive_ebsd_catalog(scanRoot);
catalog = catalog(catalog.variant == "raw", :);
[~, order] = sort(catalog.cold_reduction_percent);
catalog = catalog(order, :);
catalog = catalog([1 end], :);
assert(height(catalog) == 2);

state = repmat(empty_state(), 2, 1);
for stateIndex = 1:2
  fprintf("REFERENCE_STYLE sample=%s reduction=%.2f%%\n", ...
    catalog.sample(stateIndex), catalog.cold_reduction_percent(stateIndex));
  state(stateIndex) = calculate_state(catalog(stateIndex, :), parameters);
end

bcValues = vertcat(state.bc_values);
bcLimits = prctile(bcValues(isfinite(bcValues)), [1 99]);
grainMax = max(vertcat(state.ecd_um));
grainEdges = 0:parameters.grain_bin_width_um:...
  parameters.grain_bin_width_um * ceil(grainMax / ...
  parameters.grain_bin_width_um);
if grainEdges(end) < grainMax
  grainEdges(end + 1) = grainEdges(end) + parameters.grain_bin_width_um;
end
misorientationEdges = 0:parameters.misorientation_bin_width_deg:90;
kamEdges = 0:parameters.kam_bin_width_deg:parameters.kam_cutoff_deg;

for stateIndex = 1:2
  state(stateIndex).grain_hist_percent = weighted_histogram( ...
    state(stateIndex).ecd_um, ones(size(state(stateIndex).ecd_um)), ...
    grainEdges);
  [state(stateIndex).grain_fit_mu, state(stateIndex).grain_fit_sigma] = ...
    weighted_log_parameters(state(stateIndex).ecd_um, ...
    ones(size(state(stateIndex).ecd_um)));
  state(stateIndex).misorientation_hist_percent = weighted_histogram( ...
    state(stateIndex).boundary_angle_deg, ...
    state(stateIndex).boundary_length_um, misorientationEdges);
  state(stateIndex).kam_hist_percent = weighted_histogram( ...
    state(stateIndex).kam_values_deg, ...
    ones(size(state(stateIndex).kam_values_deg)), kamEdges);
end

pngPath = fullfile(outputDir, "reference_style_ipf_gb_kam_figure.png");
tifPath = fullfile(outputDir, "reference_style_ipf_gb_kam_figure.tif");
pdfPath = fullfile(outputDir, "reference_style_ipf_gb_kam_figure.pdf");
render_reference_figure(state, catalog, bcLimits, grainEdges, ...
  misorientationEdges, kamEdges, parameters, pngPath, tifPath, pdfPath);

summary = build_summary(state, catalog, parameters);
writetable(summary, fullfile(outputDir, ...
  "reference_style_ebsd_summary.csv"));
writetable(build_distribution_table(state, catalog, grainEdges, ...
  misorientationEdges, kamEdges), fullfile(outputDir, ...
  "reference_style_distribution_data.csv"));
writetable(struct2table(parameters), fullfile(outputDir, ...
  "reference_style_parameters.csv"));
end

function result = empty_state()
result = struct( ...
  "x", [], "y", [], "ipf_image", [], "bc_image", [], ...
  "bc_values", [], "kam_image", [], "kam_values_deg", [], ...
  "lagb_x", [], "lagb_y", [], "hagb_x", [], "hagb_y", [], ...
  "boundary_angle_deg", [], "boundary_length_um", [], ...
  "ecd_um", [], "grain_area_um2", [], "ti_crystal_symmetry", [], ...
  "grain_hist_percent", [], "grain_fit_mu", NaN, ...
  "grain_fit_sigma", NaN, "misorientation_hist_percent", [], ...
  "kam_hist_percent", [], "lagb_length_fraction", NaN, ...
  "kam_mean_deg", NaN, "ecd_mean_um", NaN, ...
  "ecd_area_weighted_mean_um", NaN, "grain_count", 0);
end

function result = calculate_state(catalogRow, parameters)
[ebsdFull, ~] = load_comprehensive_ebsd_scan(catalogRow);
tiEbsd = ebsdFull("Ti-Hex");
assert(~isempty(tiEbsd));

[xValues, yValues, nativeIndices] = native_grid_indices(ebsdFull);
result = empty_state();
result.x = xValues;
result.y = yValues;
result.ti_crystal_symmetry = tiEbsd.CS;
result.ipf_image = build_ipf_image(ebsdFull, tiEbsd, nativeIndices, ...
  numel(yValues), numel(xValues));
result.bc_values = double(ebsdFull.prop.bc(:));
result.bc_image = scalar_image(result.bc_values, nativeIndices, ...
  numel(yValues), numel(xValues));

[boundaryGrains, boundaryGrainId] = calcGrains(ebsdFull, "unitCell", ...
  "threshold", [parameters.hagb_threshold_deg ...
  parameters.boundary_detection_deg] * degree, ...
  "minPixel", parameters.min_grain_pixels);
ebsdFull.grainId = boundaryGrainId;
inner = boundaryGrains.innerBoundary("Ti-Hex", "Ti-Hex");
outer = boundaryGrains.boundary("Ti-Hex", "Ti-Hex");
innerAngle = double(angle(inner.misorientation) / degree);
outerAngle = double(angle(outer.misorientation) / degree);
lagb = inner(innerAngle >= parameters.boundary_detection_deg & ...
  innerAngle < parameters.hagb_threshold_deg);
hagb = outer(outerAngle >= parameters.hagb_threshold_deg);
[result.lagb_x, result.lagb_y] = boundary_coordinates(lagb);
[result.hagb_x, result.hagb_y] = boundary_coordinates(hagb);
result.boundary_angle_deg = [ ...
  innerAngle(innerAngle >= parameters.boundary_detection_deg & ...
  innerAngle < parameters.hagb_threshold_deg); ...
  outerAngle(outerAngle >= parameters.hagb_threshold_deg)];
result.boundary_length_um = [double(lagb.segLength(:)); ...
  double(hagb.segLength(:))];
lagbLength = sum(double(lagb.segLength));
hagbLength = sum(double(hagb.segLength));
result.lagb_length_fraction = lagbLength / (lagbLength + hagbLength);

physicalGrains = calcGrains(ebsdFull, "unitCell", ...
  "threshold", parameters.hagb_threshold_deg * degree);
tiGrains = physicalGrains("Ti-Hex");
tiGrains = tiGrains(tiGrains.numPixel >= parameters.min_grain_pixels);
grainArea = abs(double(area(tiGrains)));
grainArea = grainArea(:);
grainEcd = 2 * sqrt(grainArea / pi);
validGrain = isfinite(grainArea) & grainArea > 0 & ...
  isfinite(grainEcd) & grainEcd > 0;
result.grain_area_um2 = grainArea(validGrain);
result.ecd_um = grainEcd(validGrain);
result.grain_count = numel(result.ecd_um);
result.ecd_mean_um = mean(result.ecd_um);
result.ecd_area_weighted_mean_um = sum( ...
  result.ecd_um .* result.grain_area_um2) / sum(result.grain_area_um2);

[~, kamGrainId] = calcGrains(ebsdFull, "unitCell", ...
  "threshold", parameters.boundary_detection_deg * degree);
ebsdFull.grainId = kamGrainId;
ebsdGrid = ebsdFull.gridify;
kam = ebsdGrid.KAM("order", parameters.kam_order, ...
  "threshold", parameters.kam_cutoff_deg * degree);
kamDeg = reshape(double(kam / degree), [], 1);
tiGrid = ebsdGrid("Ti-Hex");
tiPhaseId = unique(double(tiGrid.phaseId));
assert(isscalar(tiPhaseId));
tiMask = reshape(double(ebsdGrid.phaseId), [], 1) == tiPhaseId;
kamDeg(~tiMask) = NaN;
[gridX, gridY, gridIndices] = native_grid_indices(ebsdGrid);
assert(isequal(gridX, xValues) && isequal(gridY, yValues));
result.kam_image = scalar_image(kamDeg, gridIndices, ...
  numel(gridY), numel(gridX));
result.kam_values_deg = kamDeg(isfinite(kamDeg));
result.kam_mean_deg = mean(result.kam_values_deg);
end

function render_reference_figure(state, catalog, bcLimits, grainEdges, ...
    misorientationEdges, kamEdges, parameters, pngPath, tifPath, pdfPath)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Units", "inches", "Position", [0.25 0.25 11.6 11.3]);
cleanupFigure = onCleanup(@() close(figureHandle));

mapLeft = [0.055 0.355];
mapWidth = 0.275;
mapHeight = 0.255;
rowBottom = [0.690 0.375 0.060];
statsLeft = 0.675;
statsWidth = 0.305;
statsHeight = mapHeight;
stateColors = [0.91 0.32 0.31; 0.10 0.55 0.82];

for columnIndex = 1:2
  annotation(figureHandle, "textbox", ...
    [mapLeft(columnIndex) 0.958 mapWidth 0.026], ...
    "String", sprintf("%.2f%% cold reduction | %.2f mm", ...
    catalog.cold_reduction_percent(columnIndex), ...
    catalog.diameter_mm(columnIndex)), "LineStyle", "none", ...
    "HorizontalAlignment", "center", "VerticalAlignment", "middle", ...
    "FontName", "Arial", "FontSize", 9, "FontWeight", "bold");

  ipfAxes = axes(figureHandle, "Units", "normalized", ...
    "Position", [mapLeft(columnIndex) rowBottom(1) mapWidth mapHeight]);
  imagesc(ipfAxes, state(columnIndex).x, state(columnIndex).y, ...
    state(columnIndex).ipf_image);
  format_map_axes(ipfAxes);
  draw_boundary(ipfAxes, state(columnIndex).hagb_x, ...
    state(columnIndex).hagb_y, [0.03 0.03 0.03], 0.45);
  draw_scale_bar(ipfAxes, state(columnIndex).x, state(columnIndex).y, 100);
  map_caption(ipfAxes, "IPF map");
  add_panel_label(ipfAxes, char('a' + columnIndex - 1));
  draw_ipf_inset(figureHandle, state(columnIndex).ti_crystal_symmetry, ...
    [mapLeft(columnIndex) + 0.165 rowBottom(1) + 0.165 0.095 0.075]);

  gbAxes = axes(figureHandle, "Units", "normalized", ...
    "Position", [mapLeft(columnIndex) rowBottom(2) mapWidth mapHeight]);
  imagesc(gbAxes, state(columnIndex).x, state(columnIndex).y, ...
    state(columnIndex).bc_image);
  format_map_axes(gbAxes);
  colormap(gbAxes, gray(256));
  clim(gbAxes, bcLimits);
  draw_boundary(gbAxes, state(columnIndex).lagb_x, ...
    state(columnIndex).lagb_y, [0.00 0.70 0.28], 0.50);
  draw_boundary(gbAxes, state(columnIndex).hagb_x, ...
    state(columnIndex).hagb_y, [0.03 0.03 0.03], 0.60);
  draw_scale_bar(gbAxes, state(columnIndex).x, state(columnIndex).y, 100);
  map_caption(gbAxes, "GB map");
  add_panel_label(gbAxes, char('d' + columnIndex - 1));
  draw_gb_inset_legend(gbAxes);

  kamAxes = axes(figureHandle, "Units", "normalized", ...
    "Position", [mapLeft(columnIndex) rowBottom(3) mapWidth mapHeight]);
  imageHandle = imagesc(kamAxes, state(columnIndex).x, ...
    state(columnIndex).y, state(columnIndex).kam_image);
  set(imageHandle, "AlphaData", isfinite(state(columnIndex).kam_image));
  format_map_axes(kamAxes);
  set(kamAxes, "Color", [0.78 0.78 0.78]);
  colormap(kamAxes, turbo(256));
  clim(kamAxes, [0 parameters.kam_cutoff_deg]);
  draw_boundary(kamAxes, state(columnIndex).hagb_x, ...
    state(columnIndex).hagb_y, [0.03 0.03 0.03], 0.35);
  draw_scale_bar(kamAxes, state(columnIndex).x, state(columnIndex).y, 100);
  map_caption(kamAxes, "KAM map");
  add_panel_label(kamAxes, char('g' + columnIndex - 1));
  draw_horizontal_color_scale(figureHandle, ...
    [mapLeft(columnIndex) + 0.045 rowBottom(3) + mapHeight + 0.008 ...
    mapWidth - 0.09 0.015], [0 parameters.kam_cutoff_deg]);
end

draw_paired_grain_histograms(figureHandle, ...
  [statsLeft rowBottom(1) statsWidth statsHeight], state, catalog, ...
  grainEdges, stateColors);
draw_paired_misorientation_histograms(figureHandle, ...
  [statsLeft rowBottom(2) statsWidth statsHeight], state, catalog, ...
  misorientationEdges, stateColors, parameters.hagb_threshold_deg);
draw_paired_kam_histograms(figureHandle, ...
  [statsLeft rowBottom(3) statsWidth statsHeight], state, catalog, ...
  kamEdges, stateColors);

annotation(figureHandle, "textbox", [statsLeft - 0.020 0.930 0.035 0.025], ...
  "String", "(c)", "LineStyle", "none", "FontName", "Arial", ...
  "FontSize", 12, "FontWeight", "bold");
annotation(figureHandle, "textbox", [statsLeft - 0.020 0.615 0.035 0.025], ...
  "String", "(f)", "LineStyle", "none", "FontName", "Arial", ...
  "FontSize", 12, "FontWeight", "bold");
annotation(figureHandle, "textbox", [statsLeft - 0.020 0.300 0.035 0.025], ...
  "String", "(i)", "LineStyle", "none", "FontName", "Arial", ...
  "FontSize", 12, "FontWeight", "bold");

drawnow;
exportgraphics(figureHandle, pngPath, "Resolution", 600, ...
  "BackgroundColor", "white");
exportgraphics(figureHandle, tifPath, "Resolution", 600, ...
  "BackgroundColor", "white");
exportgraphics(figureHandle, pdfPath, "ContentType", "image", ...
  "Resolution", 600, "BackgroundColor", "white");
clear cleanupFigure
end

function draw_paired_grain_histograms(figureHandle, position, state, ...
    catalog, edges, colors)
axesHandles = paired_axes(figureHandle, position);
centers = edges(1:end-1) + diff(edges) / 2;
fineX = linspace(max(0.05, edges(1) + 0.05), edges(end), 500);
binWidth = median(diff(edges));
for stateIndex = 1:2
  axesHandle = axesHandles(stateIndex);
  bar(axesHandle, centers, state(stateIndex).grain_hist_percent, 0.82, ...
    "FaceColor", colors(stateIndex, :), "FaceAlpha", 0.58, ...
    "EdgeColor", [0.25 0.25 0.25], "LineWidth", 0.35);
  hold(axesHandle, "on");
  fitValues = 100 * binWidth * lognormal_pdf(fineX, ...
    state(stateIndex).grain_fit_mu, state(stateIndex).grain_fit_sigma);
  plot(axesHandle, fineX, fitValues, "Color", colors(stateIndex, :), ...
    "LineWidth", 1.6);
  yMax = max([state(stateIndex).grain_hist_percent(:); fitValues(:)]);
  ylim(axesHandle, [0 max(5, ceil(1.12 * yMax / 2) * 2)]);
  xlim(axesHandle, [edges(1) edges(end)]);
  format_distribution_axes(axesHandle, stateIndex == 2);
  text(axesHandle, 0.02, 0.86, sprintf("Avg. = %.2f um", ...
    state(stateIndex).ecd_mean_um), "Units", "normalized", ...
    "Color", colors(stateIndex, :), "FontName", "Arial", ...
    "FontSize", 8.5, "FontWeight", "bold", "VerticalAlignment", "top");
  state_tag(axesHandle, catalog(stateIndex, :), colors(stateIndex, :));
end
xlabel(axesHandles(2), "Grain size (um)", "FontName", "Arial", ...
  "FontSize", 8);
ylabel(axesHandles(1), "Frequency (%)", "FontName", "Arial", ...
  "FontSize", 8);
ylabel(axesHandles(2), "Frequency (%)", "FontName", "Arial", ...
  "FontSize", 8);
end

function draw_paired_misorientation_histograms(figureHandle, position, ...
    state, catalog, edges, colors, hagbThresholdDeg)
axesHandles = paired_axes(figureHandle, position);
centers = edges(1:end-1) + diff(edges) / 2;
for stateIndex = 1:2
  axesHandle = axesHandles(stateIndex);
  bar(axesHandle, centers, ...
    state(stateIndex).misorientation_hist_percent, 0.82, ...
    "FaceColor", colors(stateIndex, :), "FaceAlpha", 0.58, ...
    "EdgeColor", [0.25 0.25 0.25], "LineWidth", 0.30);
  hold(axesHandle, "on");
  xline(axesHandle, hagbThresholdDeg, "--", "Color", [0.1 0.1 0.1], ...
    "LineWidth", 1.0);
  yMax = max(state(stateIndex).misorientation_hist_percent);
  ylim(axesHandle, [0 max(5, ceil(1.22 * yMax / 5) * 5)]);
  xlim(axesHandle, [edges(1) edges(end)]);
  format_distribution_axes(axesHandle, stateIndex == 2);
  text(axesHandle, 0.98, 0.84, sprintf("LAGBs = %.1f%%", ...
    100 * state(stateIndex).lagb_length_fraction), ...
    "Units", "normalized", "HorizontalAlignment", "right", ...
    "VerticalAlignment", "top", "Color", colors(stateIndex, :), ...
    "FontName", "Arial", "FontSize", 8.5, "FontWeight", "bold");
  state_tag(axesHandle, catalog(stateIndex, :), colors(stateIndex, :));
end
xlabel(axesHandles(2), "Misorientation angle (deg)", ...
  "FontName", "Arial", "FontSize", 8);
ylabel(axesHandles(1), "Frequency (%)", "FontName", "Arial", ...
  "FontSize", 8);
ylabel(axesHandles(2), "Frequency (%)", "FontName", "Arial", ...
  "FontSize", 8);
end

function draw_paired_kam_histograms(figureHandle, position, state, ...
    catalog, edges, colors)
axesHandles = paired_axes(figureHandle, position);
centers = edges(1:end-1) + diff(edges) / 2;
for stateIndex = 1:2
  axesHandle = axesHandles(stateIndex);
  bar(axesHandle, centers, state(stateIndex).kam_hist_percent, 0.82, ...
    "FaceColor", colors(stateIndex, :), "FaceAlpha", 0.62, ...
    "EdgeColor", [0.25 0.25 0.25], "LineWidth", 0.30);
  yMaximum = max(5, ceil(1.18 * ...
    max(state(stateIndex).kam_hist_percent) / 5) * 5);
  ylim(axesHandle, [0 yMaximum]);
  xlim(axesHandle, [edges(1) edges(end)]);
  format_distribution_axes(axesHandle, stateIndex == 2);
  text(axesHandle, 0.98, 0.84, sprintf("KAM = %.2f deg", ...
    state(stateIndex).kam_mean_deg), "Units", "normalized", ...
    "HorizontalAlignment", "right", "VerticalAlignment", "top", ...
    "Color", colors(stateIndex, :), "FontName", "Arial", ...
    "FontSize", 8.5, "FontWeight", "bold");
  state_tag(axesHandle, catalog(stateIndex, :), colors(stateIndex, :));
end
xlabel(axesHandles(2), "Kernel average misorientation (deg)", ...
  "FontName", "Arial", "FontSize", 8);
ylabel(axesHandles(1), "Frequency (%)", "FontName", "Arial", ...
  "FontSize", 8);
ylabel(axesHandles(2), "Frequency (%)", "FontName", "Arial", ...
  "FontSize", 8);
end

function axesHandles = paired_axes(figureHandle, position)
gap = 0.006;
subHeight = (position(4) - gap) / 2;
axesHandles = gobjects(2, 1);
axesHandles(1) = axes(figureHandle, "Units", "normalized", ...
  "Position", [position(1) position(2) + subHeight + gap ...
  position(3) subHeight]);
axesHandles(2) = axes(figureHandle, "Units", "normalized", ...
  "Position", [position(1) position(2) position(3) subHeight]);
end

function format_distribution_axes(axesHandle, showXLabels)
set(axesHandle, "FontName", "Arial", "FontSize", 7, ...
  "TickDir", "out", "Box", "on", "LineWidth", 0.6, "Layer", "top");
if ~showXLabels
  set(axesHandle, "XTickLabel", []);
end
end

function state_tag(axesHandle, catalogRow, color)
text(axesHandle, 0.98, 0.98, sprintf("%.2f%%", ...
  catalogRow.cold_reduction_percent), "Units", "normalized", ...
  "HorizontalAlignment", "right", "VerticalAlignment", "top", ...
  "FontName", "Arial", "FontSize", 6.8, "Color", color);
end

function draw_ipf_inset(figureHandle, crystalSymmetry, position)
axesHandle = axes(figureHandle, "Units", "normalized", ...
  "Position", position, "Color", "white");
colorKey = ipfHSVKey(crystalSymmetry);
colorKey.inversePoleFigureDirection = xvector;
plot(colorKey, "parent", axesHandle, "noTitle");
end

function draw_gb_inset_legend(axesHandle)
xLimits = xlim(axesHandle);
yLimits = ylim(axesHandle);
xRange = diff(xLimits);
yRange = diff(yLimits);
line(axesHandle, xLimits(1) + [0.05 0.17] * xRange, ...
  yLimits(1) + [0.94 0.94] * yRange, ...
  "Color", [0.03 0.03 0.03], "LineWidth", 1.5);
text(axesHandle, 0.19, 0.94, "HAGBs", "Units", "normalized", ...
  "VerticalAlignment", "middle", "FontName", "Arial", ...
  "FontSize", 7, "FontWeight", "bold", "BackgroundColor", "white");
line(axesHandle, xLimits(1) + [0.50 0.62] * xRange, ...
  yLimits(1) + [0.94 0.94] * yRange, ...
  "Color", [0.00 0.70 0.28], "LineWidth", 1.5);
text(axesHandle, 0.64, 0.94, "LAGBs", "Units", "normalized", ...
  "VerticalAlignment", "middle", "FontName", "Arial", ...
  "FontSize", 7, "FontWeight", "bold", "BackgroundColor", "white", ...
  "Color", [0.00 0.55 0.20]);
end

function draw_horizontal_color_scale(figureHandle, position, limits)
axesHandle = axes(figureHandle, "Units", "normalized", ...
  "Position", position);
imagesc(axesHandle, linspace(limits(1), limits(2), 256), 1, 1:256);
set(axesHandle, "YTick", [], "XTick", limits, "XLim", limits, ...
  "FontName", "Arial", "FontSize", 6, "Box", "on", ...
  "XAxisLocation", "top");
colormap(axesHandle, turbo(256));
end

function map_caption(axesHandle, caption)
text(axesHandle, 0.98, 0.035, caption, "Units", "normalized", ...
  "HorizontalAlignment", "right", "VerticalAlignment", "bottom", ...
  "FontName", "Arial", "FontSize", 9, "FontWeight", "bold", ...
  "BackgroundColor", "white", "Margin", 1);
end

function add_panel_label(axesHandle, letter)
text(axesHandle, 0.012, 0.985, "(" + string(letter) + ")", ...
  "Units", "normalized", "HorizontalAlignment", "left", ...
  "VerticalAlignment", "top", "FontName", "Arial", ...
  "FontSize", 12, "FontWeight", "bold", "BackgroundColor", "white", ...
  "Margin", 1);
end

function format_map_axes(axesHandle)
axis(axesHandle, "image");
set(axesHandle, "YDir", "normal", "XTick", [], "YTick", [], ...
  "Box", "on", "LineWidth", 0.65, "FontName", "Arial");
hold(axesHandle, "on");
end

function draw_boundary(axesHandle, xCoordinates, yCoordinates, color, width)
line(axesHandle, xCoordinates, yCoordinates, "Color", color, ...
  "LineWidth", width);
end

function draw_scale_bar(axesHandle, xValues, yValues, lengthUm)
xRange = max(xValues) - min(xValues);
yRange = max(yValues) - min(yValues);
xStart = min(xValues) + 0.055 * xRange;
yPosition = min(yValues) + 0.065 * yRange;
line(axesHandle, [xStart xStart + lengthUm], [yPosition yPosition], ...
  "Color", "black", "LineWidth", 2.1);
text(axesHandle, xStart + lengthUm / 2, yPosition + 0.018 * yRange, ...
  sprintf("%d um", lengthUm), "HorizontalAlignment", "center", ...
  "VerticalAlignment", "bottom", "FontName", "Arial", ...
  "FontSize", 7, "BackgroundColor", "white", "Margin", 0.6);
end

function imageData = build_ipf_image(ebsdFull, tiEbsd, nativeIndices, ...
    rowCount, columnCount)
imageData = repmat(reshape([0.82 0.82 0.82], 1, 1, 3), ...
  rowCount, columnCount, 1);
colorKey = ipfHSVKey(tiEbsd);
colorKey.inversePoleFigureDirection = xvector;
colors = colorKey.orientation2color(tiEbsd.orientations);
[found, fullRows] = ismember(double(tiEbsd.id(:)), double(ebsdFull.id(:)));
assert(all(found));
for channelIndex = 1:3
  channel = imageData(:, :, channelIndex);
  channel(nativeIndices(fullRows)) = colors(:, channelIndex);
  imageData(:, :, channelIndex) = channel;
end
end

function imageData = scalar_image(values, indices, rowCount, columnCount)
values = reshape(double(values), [], 1);
assert(numel(values) == numel(indices));
imageData = nan(rowCount, columnCount);
imageData(indices) = values;
end

function [xValues, yValues, linearIndices] = native_grid_indices(ebsd)
xCoordinates = double(ebsd.x(:));
yCoordinates = double(ebsd.y(:));
xValues = unique(xCoordinates);
yValues = unique(yCoordinates);
[xFound, xIndex] = ismember(xCoordinates, xValues);
[yFound, yIndex] = ismember(yCoordinates, yValues);
assert(all(xFound & yFound));
linearIndices = sub2ind([numel(yValues), numel(xValues)], yIndex, xIndex);
assert(numel(unique(linearIndices)) == length(ebsd));
end

function [xCoordinates, yCoordinates] = boundary_coordinates(boundary)
if isempty(boundary)
  xCoordinates = nan;
  yCoordinates = nan;
  return
end
vertices = boundary.allV.xyz;
faces = boundary.F;
xCoordinates = reshape(vertices(faces.', 1), 2, []);
yCoordinates = reshape(vertices(faces.', 2), 2, []);
separator = nan(1, size(xCoordinates, 2));
xCoordinates = reshape([xCoordinates; separator], [], 1);
yCoordinates = reshape([yCoordinates; separator], [], 1);
end

function percentages = weighted_histogram(values, weights, edges)
values = double(values(:));
weights = double(weights(:));
valid = isfinite(values) & isfinite(weights) & weights > 0 & ...
  values >= edges(1) & values <= edges(end);
values = values(valid);
weights = weights(valid);
bins = discretize(values, edges);
assert(all(isfinite(bins)));
percentages = 100 * accumarray(bins, weights, ...
  [numel(edges) - 1 1], @sum, 0).' / sum(weights);
end

function [mu, sigma] = weighted_log_parameters(values, weights)
values = double(values(:));
weights = double(weights(:));
valid = isfinite(values) & values > 0 & isfinite(weights) & weights > 0;
values = values(valid);
weights = weights(valid) / sum(weights(valid));
logValues = log(values);
mu = sum(weights .* logValues);
sigma = sqrt(sum(weights .* (logValues - mu).^2));
sigma = max(sigma, sqrt(eps));
end

function density = lognormal_pdf(values, mu, sigma)
density = exp(-0.5 * ((log(values) - mu) / sigma).^2) ./ ...
  (values * sigma * sqrt(2 * pi));
end

function summary = build_summary(state, catalog, parameters)
sample = catalog.sample;
diameter_mm = catalog.diameter_mm;
cold_reduction_percent = catalog.cold_reduction_percent;
variant = catalog.variant;
grain_count = reshape([state.grain_count], [], 1);
ecd_number_mean_um = reshape([state.ecd_mean_um], [], 1);
ecd_area_weighted_mean_um = reshape( ...
  [state.ecd_area_weighted_mean_um], [], 1);
lagb_length_fraction = reshape([state.lagb_length_fraction], [], 1);
kam_mean_deg = reshape([state.kam_mean_deg], [], 1);
boundary_detection_deg = repmat(parameters.boundary_detection_deg, 2, 1);
hagb_threshold_deg = repmat(parameters.hagb_threshold_deg, 2, 1);
min_grain_pixels = repmat(parameters.min_grain_pixels, 2, 1);
kam_order = parameters.kam_order * ones(2, 1);
kam_cutoff_deg = repmat(parameters.kam_cutoff_deg, 2, 1);
summary = table(sample, diameter_mm, cold_reduction_percent, variant, ...
  grain_count, ecd_number_mean_um, ecd_area_weighted_mean_um, ...
  lagb_length_fraction, kam_mean_deg, boundary_detection_deg, ...
  hagb_threshold_deg, min_grain_pixels, kam_order, kam_cutoff_deg);
end

function output = build_distribution_table(state, catalog, grainEdges, ...
    misorientationEdges, kamEdges)
plot_type = strings(0, 1);
sample = strings(0, 1);
cold_reduction_percent = zeros(0, 1);
bin_lower = zeros(0, 1);
bin_upper = zeros(0, 1);
bin_center = zeros(0, 1);
frequency_percent = zeros(0, 1);
for stateIndex = 1:2
  [plot_type, sample, cold_reduction_percent, bin_lower, bin_upper, ...
    bin_center, frequency_percent] = append_distribution( ...
    plot_type, sample, cold_reduction_percent, bin_lower, bin_upper, ...
    bin_center, frequency_percent, "grain_size_number", ...
    catalog(stateIndex, :), grainEdges, ...
    state(stateIndex).grain_hist_percent);
  [plot_type, sample, cold_reduction_percent, bin_lower, bin_upper, ...
    bin_center, frequency_percent] = append_distribution( ...
    plot_type, sample, cold_reduction_percent, bin_lower, bin_upper, ...
    bin_center, frequency_percent, "boundary_length_weighted", ...
    catalog(stateIndex, :), misorientationEdges, ...
    state(stateIndex).misorientation_hist_percent);
  [plot_type, sample, cold_reduction_percent, bin_lower, bin_upper, ...
    bin_center, frequency_percent] = append_distribution( ...
    plot_type, sample, cold_reduction_percent, bin_lower, bin_upper, ...
    bin_center, frequency_percent, "kam_pixel_frequency", ...
    catalog(stateIndex, :), kamEdges, state(stateIndex).kam_hist_percent);
end
output = table(plot_type, sample, cold_reduction_percent, bin_lower, ...
  bin_upper, bin_center, frequency_percent);
end

function [plotTypeOut, sampleOut, reductionOut, lowerOut, upperOut, ...
    centerOut, frequencyOut] = append_distribution(plotTypeIn, sampleIn, ...
    reductionIn, lowerIn, upperIn, centerIn, frequencyIn, plotType, ...
    catalogRow, edges, frequency)
n = numel(frequency);
plotTypeOut = [plotTypeIn; repmat(plotType, n, 1)];
sampleOut = [sampleIn; repmat(catalogRow.sample, n, 1)];
reductionOut = [reductionIn; repmat( ...
  catalogRow.cold_reduction_percent, n, 1)];
lowerOut = [lowerIn; edges(1:end-1).'];
upperOut = [upperIn; edges(2:end).'];
centerOut = [centerIn; (edges(1:end-1) + diff(edges) / 2).'];
frequencyOut = [frequencyIn; frequency(:)];
end
