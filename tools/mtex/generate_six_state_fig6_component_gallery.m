function summary = generate_six_state_fig6_component_gallery(scanRoot, outputDir)
%GENERATE_SIX_STATE_FIG6_COMPONENT_GALLERY Export six EBSD component sets.
% Each component is exported as six standalone 600-dpi PNG panels and one
% 2-by-3 comparison montage. Grain-size panels report 2 degree
% orientation-domain equivalent-circle diameters (ECD), not 15 degree
% HAGB-defined physical grain sizes.

arguments
  scanRoot (1,1) string
  outputDir (1,1) string
end

assert(isfolder(scanRoot), "EBSD scan folder not found: %s", scanRoot);
assert(~isempty(which("EBSD")), "MTEX must be loaded before this function.");

parameters = struct( ...
  "boundary_detection_deg", 2, ...
  "hagb_threshold_deg", 15, ...
  "min_domain_pixels", 5, ...
  "kam_order", 1, ...
  "kam_cutoff_deg", 5, ...
  "kam_display_max_deg", 5, ...
  "domain_size_bin_width_um", 2, ...
  "boundary_bin_width_deg", 2, ...
  "kam_bin_width_deg", 0.2, ...
  "scale_bar_um", 100, ...
  "export_dpi", 600);

folders = create_output_folders(outputDir);
catalog = comprehensive_ebsd_catalog(scanRoot);
catalog = catalog(catalog.variant == "raw", :);
[~, order] = sort(catalog.cold_reduction_percent);
catalog = catalog(order, :);
assert(height(catalog) == 6, "Expected six registered raw EBSD states.");

state = repmat(empty_state(), height(catalog), 1);
for stateIndex = 1:height(catalog)
  fprintf("COMPONENT_GALLERY sample=%s reduction=%.2f%%\n", ...
    catalog.sample(stateIndex), catalog.cold_reduction_percent(stateIndex));
  state(stateIndex) = calculate_state(catalog(stateIndex, :), parameters);
end

allBc = vertcat(state.bc_values);
bcLimits = robust_limits(allBc);
allEcd = vertcat(state.domain_ecd_um);
domainEdges = common_edges(allEcd, parameters.domain_size_bin_width_um, 0);
allBoundaryAngles = vertcat(state.boundary_angle_deg);
boundaryEdges = common_edges(allBoundaryAngles, ...
  parameters.boundary_bin_width_deg, 0, 90);
kamEdges = 0:parameters.kam_bin_width_deg:parameters.kam_cutoff_deg;

for stateIndex = 1:height(catalog)
  state(stateIndex).domain_hist_percent = weighted_histogram( ...
    state(stateIndex).domain_ecd_um, ...
    ones(size(state(stateIndex).domain_ecd_um)), domainEdges);
  [state(stateIndex).domain_fit_mu, state(stateIndex).domain_fit_sigma] = ...
    weighted_log_parameters(state(stateIndex).domain_ecd_um, ...
    ones(size(state(stateIndex).domain_ecd_um)));
  state(stateIndex).boundary_hist_percent = weighted_histogram( ...
    state(stateIndex).boundary_angle_deg, ...
    state(stateIndex).boundary_length_um, boundaryEdges);
  state(stateIndex).kam_hist_percent = weighted_histogram( ...
    state(stateIndex).kam_values_deg, ...
    ones(size(state(stateIndex).kam_values_deg)), kamEdges);
end

colors = [ ...
  0.40 0.40 0.40
  0.00 0.45 0.70
  0.34 0.71 0.91
  0.00 0.62 0.45
  0.90 0.62 0.00
  0.84 0.37 0.00];

export_standalone_maps(state, catalog, folders, bcLimits, parameters);
export_map_montages(state, catalog, folders, bcLimits, parameters);
export_standalone_distributions(state, catalog, folders, domainEdges, ...
  boundaryEdges, kamEdges, colors, parameters);
export_distribution_montages(state, catalog, folders, domainEdges, ...
  boundaryEdges, kamEdges, colors, parameters);

summary = build_summary(state, catalog, parameters);
writetable(summary, fullfile(outputDir, "component_summary.csv"));
writetable(build_distribution_table(state, catalog, domainEdges, ...
  boundaryEdges, kamEdges), fullfile(outputDir, "distribution_data.csv"));
writetable(struct2table(parameters), fullfile(outputDir, "parameters.csv"));

validate_outputs(state, domainEdges, boundaryEdges, kamEdges, folders);
end

function folders = create_output_folders(outputDir)
folders = struct( ...
  "ipf", fullfile(outputDir, "01_ipf_maps"), ...
  "gb", fullfile(outputDir, "02_gb_maps"), ...
  "kam", fullfile(outputDir, "03_kam_maps"), ...
  "domain", fullfile(outputDir, "04_orientation_domain_size"), ...
  "boundary", fullfile(outputDir, "05_boundary_misorientation"), ...
  "kamdist", fullfile(outputDir, "06_kam_distributions"), ...
  "montage", fullfile(outputDir, "montages"));
names = fieldnames(folders);
for index = 1:numel(names)
  if ~isfolder(folders.(names{index}))
    mkdir(folders.(names{index}));
  end
end
end

function result = empty_state()
result = struct( ...
  "x", [], "y", [], "ipf_image", [], "bc_image", [], ...
  "bc_values", [], "kam_image", [], "kam_values_deg", [], ...
  "lagb_x", [], "lagb_y", [], "hagb_x", [], "hagb_y", [], ...
  "boundary_angle_deg", [], "boundary_length_um", [], ...
  "domain_ecd_um", [], "domain_area_um2", [], ...
  "ti_crystal_symmetry", [], "domain_hist_percent", [], ...
  "domain_fit_mu", NaN, "domain_fit_sigma", NaN, ...
  "boundary_hist_percent", [], "kam_hist_percent", [], ...
  "domain_count", 0, "domain_mean_um", NaN, ...
  "domain_median_um", NaN, "domain_area_weighted_mean_um", NaN, ...
  "domain_area_weighted_median_um", NaN, ...
  "lagb_length_fraction", NaN, "hagb_length_fraction", NaN, ...
  "kam_mean_deg", NaN, "kam_median_deg", NaN, "kam_p90_deg", NaN);
end

function result = calculate_state(catalogRow, parameters)
[ebsdFull, ~] = load_comprehensive_ebsd_scan(catalogRow);
tiEbsd = ebsdFull("Ti-Hex");
assert(~isempty(tiEbsd), "Ti-Hex phase is absent.");
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
  "minPixel", parameters.min_domain_pixels);
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
totalBoundaryLength = lagbLength + hagbLength;
result.lagb_length_fraction = lagbLength / totalBoundaryLength;
result.hagb_length_fraction = hagbLength / totalBoundaryLength;

domains = calcGrains(ebsdFull, "unitCell", ...
  "threshold", parameters.boundary_detection_deg * degree);
domains = domains("Ti-Hex");
domains = domains(domains.numPixel >= parameters.min_domain_pixels);
domainArea = abs(double(area(domains)));
domainArea = domainArea(:);
domainEcd = 2 * sqrt(domainArea / pi);
validDomain = isfinite(domainArea) & domainArea > 0 & ...
  isfinite(domainEcd) & domainEcd > 0;
result.domain_area_um2 = domainArea(validDomain);
result.domain_ecd_um = domainEcd(validDomain);
assert(~isempty(result.domain_ecd_um), "No 2 degree domains pass minPixel.");
result.domain_count = numel(result.domain_ecd_um);
result.domain_mean_um = mean(result.domain_ecd_um);
result.domain_median_um = median(result.domain_ecd_um);
result.domain_area_weighted_mean_um = sum(result.domain_ecd_um .* ...
  result.domain_area_um2) / sum(result.domain_area_um2);
result.domain_area_weighted_median_um = weighted_quantile( ...
  result.domain_ecd_um, result.domain_area_um2, 0.5);

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
assert(isequal(gridX, xValues) && isequal(gridY, yValues), ...
  "Gridify changed registered map coordinates.");
result.kam_image = scalar_image(kamDeg, gridIndices, ...
  numel(gridY), numel(gridX));
result.kam_values_deg = kamDeg(isfinite(kamDeg));
result.kam_mean_deg = mean(result.kam_values_deg);
result.kam_median_deg = median(result.kam_values_deg);
result.kam_p90_deg = prctile(result.kam_values_deg, 90);
end

function export_standalone_maps(state, catalog, folders, bcLimits, parameters)
for stateIndex = 1:numel(state)
  tag = state_tag(catalog(stateIndex, :));
  export_one_map(state(stateIndex), catalog(stateIndex, :), "ipf", ...
    bcLimits, parameters, fullfile(folders.ipf, tag + "_ipf_ad.png"));
  export_one_map(state(stateIndex), catalog(stateIndex, :), "gb", ...
    bcLimits, parameters, fullfile(folders.gb, tag + "_gb.png"));
  export_one_map(state(stateIndex), catalog(stateIndex, :), "kam", ...
    bcLimits, parameters, fullfile(folders.kam, tag + "_kam.png"));
end
end

function export_one_map(state, catalogRow, mapType, bcLimits, parameters, path)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Units", "inches", "Position", [0.25 0.25 4.2 4.6]);
cleanupFigure = onCleanup(@() close(figureHandle));
switch mapType
  case "ipf"
    axesPosition = [0.075 0.205 0.85 0.675];
  case "kam"
    axesPosition = [0.075 0.060 0.85 0.735];
  otherwise
    axesPosition = [0.075 0.075 0.85 0.82];
end
axesHandle = axes(figureHandle, "Units", "normalized", ...
  "Position", axesPosition);
draw_map_panel(axesHandle, state, mapType, bcLimits, parameters);
title(axesHandle, sprintf("%.2f%% cold reduction | %.2f mm", ...
  catalogRow.cold_reduction_percent, catalogRow.diameter_mm), ...
  "FontName", "Arial", "FontSize", 9, "FontWeight", "bold");
switch mapType
  case "ipf"
    draw_ipf_inset(figureHandle, state.ti_crystal_symmetry, ...
      [0.38 0.075 0.24 0.090]);
  case "gb"
    draw_gb_inset_legend(axesHandle);
  case "kam"
    draw_horizontal_color_scale(figureHandle, [0.23 0.905 0.54 0.022], ...
      [0 parameters.kam_display_max_deg]);
end
exportgraphics(figureHandle, path, "Resolution", parameters.export_dpi, ...
  "BackgroundColor", "white");
clear cleanupFigure
end

function export_map_montages(state, catalog, folders, bcLimits, parameters)
types = ["ipf", "gb", "kam"];
names = ["ipf_ad_six_state_montage", "gb_six_state_montage", ...
  "kam_six_state_montage"];
for typeIndex = 1:numel(types)
  figureHandle = figure("Visible", "off", "Color", "white", ...
    "Units", "inches", "Position", [0.25 0.25 11.2 7.5]);
  cleanupFigure = onCleanup(@() close(figureHandle));
  left = [0.035 0.355 0.675];
  bottom = [0.525 0.105];
  panelWidth = 0.29;
  panelHeight = 0.355;
  for stateIndex = 1:6
    rowIndex = 1 + floor((stateIndex - 1) / 3);
    columnIndex = 1 + mod(stateIndex - 1, 3);
    axesHandle = axes(figureHandle, "Units", "normalized", ...
      "Position", [left(columnIndex) bottom(rowIndex) panelWidth panelHeight]);
    draw_map_panel(axesHandle, state(stateIndex), types(typeIndex), ...
      bcLimits, parameters);
    title(axesHandle, sprintf("(%c) %.2f%% | %.2f mm", ...
      char('a' + stateIndex - 1), ...
      catalog.cold_reduction_percent(stateIndex), ...
      catalog.diameter_mm(stateIndex)), "FontName", "Arial", ...
      "FontSize", 8, "FontWeight", "bold");
  end
  switch types(typeIndex)
    case "ipf"
      draw_ipf_inset(figureHandle, state(1).ti_crystal_symmetry, ...
        [0.835 0.005 0.11 0.08]);
    case "gb"
      draw_montage_gb_legend(figureHandle, [0.34 0.012 0.32 0.055]);
    case "kam"
      draw_horizontal_color_scale(figureHandle, [0.36 0.025 0.28 0.025], ...
        [0 parameters.kam_display_max_deg]);
  end
  annotation(figureHandle, "textbox", [0.01 0.965 0.98 0.025], ...
    "String", montage_title(types(typeIndex)), "LineStyle", "none", ...
    "HorizontalAlignment", "center", "FontName", "Arial", ...
    "FontSize", 10, "FontWeight", "bold");
  export_figure_triplet(figureHandle, folders.montage, names(typeIndex), ...
    parameters.export_dpi, "image");
  clear cleanupFigure
end
end

function value = montage_title(mapType)
switch mapType
  case "ipf"
    value = "IPF-AD maps | raw EBSD";
  case "gb"
    value = "Grain-boundary maps | 2-<15 deg LAGB; >=15 deg HAGB";
  otherwise
    value = "KAM maps | first neighbour; 5 deg cutoff; common 0-5 deg scale";
end
end

function draw_map_panel(axesHandle, state, mapType, bcLimits, parameters)
switch mapType
  case "ipf"
    imagesc(axesHandle, state.x, state.y, state.ipf_image);
    format_map_axes(axesHandle);
    draw_boundary(axesHandle, state.hagb_x, state.hagb_y, ...
      [0.03 0.03 0.03], 0.38);
  case "gb"
    imagesc(axesHandle, state.x, state.y, state.bc_image);
    format_map_axes(axesHandle);
    colormap(axesHandle, gray(256));
    clim(axesHandle, bcLimits);
    draw_boundary(axesHandle, state.lagb_x, state.lagb_y, ...
      [0.00 0.70 0.28], 0.55);
    draw_boundary(axesHandle, state.hagb_x, state.hagb_y, ...
      [0.03 0.03 0.03], 0.65);
  case "kam"
    imageHandle = imagesc(axesHandle, state.x, state.y, state.kam_image);
    set(imageHandle, "AlphaData", isfinite(state.kam_image));
    format_map_axes(axesHandle);
    set(axesHandle, "Color", [0.80 0.80 0.80]);
    colormap(axesHandle, turbo(256));
    clim(axesHandle, [0 parameters.kam_display_max_deg]);
    draw_boundary(axesHandle, state.hagb_x, state.hagb_y, ...
      [0.03 0.03 0.03], 0.32);
    text(axesHandle, 0.98, 0.04, sprintf("Mean = %.3f deg", ...
      state.kam_mean_deg), "Units", "normalized", ...
      "HorizontalAlignment", "right", "VerticalAlignment", "bottom", ...
      "FontName", "Arial", "FontSize", 7, "BackgroundColor", "white", ...
      "Margin", 1);
end
draw_scale_bar(axesHandle, state.x, state.y, parameters.scale_bar_um);
end

function export_standalone_distributions(state, catalog, folders, ...
    domainEdges, boundaryEdges, kamEdges, colors, parameters)
for stateIndex = 1:6
  tag = state_tag(catalog(stateIndex, :));
  domainYMax = domain_ymax(state(stateIndex), domainEdges);
  boundaryYMax = common_ymax(state(stateIndex).boundary_hist_percent);
  kamYMax = common_ymax(state(stateIndex).kam_hist_percent);
  export_one_distribution(state(stateIndex), catalog(stateIndex, :), ...
    "domain", domainEdges, colors(stateIndex, :), domainYMax, ...
    parameters, fullfile(folders.domain, tag + "_orientation_domain_ecd.png"));
  export_one_distribution(state(stateIndex), catalog(stateIndex, :), ...
    "boundary", boundaryEdges, colors(stateIndex, :), boundaryYMax, ...
    parameters, fullfile(folders.boundary, tag + "_boundary_misorientation.png"));
  export_one_distribution(state(stateIndex), catalog(stateIndex, :), ...
    "kam", kamEdges, colors(stateIndex, :), kamYMax, parameters, ...
    fullfile(folders.kamdist, tag + "_kam_distribution.png"));
end
end

function export_one_distribution(state, catalogRow, plotType, edges, ...
    color, yMaximum, parameters, path)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Units", "inches", "Position", [0.25 0.25 4.2 3.2]);
cleanupFigure = onCleanup(@() close(figureHandle));
  axesHandle = axes(figureHandle, "Units", "normalized", ...
  "Position", [0.17 0.18 0.78 0.69]);
draw_distribution_panel(axesHandle, state, plotType, edges, color, ...
  yMaximum, parameters, true, true);
title(axesHandle, sprintf("%.2f%% cold reduction | %.2f mm", ...
  catalogRow.cold_reduction_percent, catalogRow.diameter_mm), ...
  "FontName", "Arial", "FontSize", 9, "FontWeight", "bold");
exportgraphics(figureHandle, path, "Resolution", parameters.export_dpi, ...
  "BackgroundColor", "white");
clear cleanupFigure
end

function export_distribution_montages(state, catalog, folders, ...
    domainEdges, boundaryEdges, kamEdges, colors, parameters)
types = ["domain", "boundary", "kam"];
names = ["orientation_domain_ecd_six_state_montage", ...
  "boundary_misorientation_six_state_montage", ...
  "kam_distribution_six_state_montage"];
edgeSets = {domainEdges, boundaryEdges, kamEdges};
frequencySets = {vertcat(state.domain_hist_percent), ...
  vertcat(state.boundary_hist_percent), vertcat(state.kam_hist_percent)};
for typeIndex = 1:3
  figureHandle = figure("Visible", "off", "Color", "white", ...
    "Units", "inches", "Position", [0.25 0.25 11.2 7.2]);
  cleanupFigure = onCleanup(@() close(figureHandle));
  layout = tiledlayout(figureHandle, 2, 3, "Padding", "compact", ...
    "TileSpacing", "compact");
  if types(typeIndex) == "domain"
    yMaximum = domain_ymax(state, edgeSets{typeIndex});
  else
    yMaximum = common_ymax(frequencySets{typeIndex});
  end
  for stateIndex = 1:6
    axesHandle = nexttile(layout);
    draw_distribution_panel(axesHandle, state(stateIndex), ...
      types(typeIndex), edgeSets{typeIndex}, colors(stateIndex, :), ...
      yMaximum, parameters, true, stateIndex == 1);
    title(axesHandle, sprintf("(%c) %.2f%% | %.2f mm", ...
      char('a' + stateIndex - 1), ...
      catalog.cold_reduction_percent(stateIndex), ...
      catalog.diameter_mm(stateIndex)), "FontName", "Arial", ...
      "FontSize", 8, "FontWeight", "bold");
  end
  title(layout, distribution_montage_title(types(typeIndex)), ...
    "FontName", "Arial", "FontSize", 10, "FontWeight", "bold");
  export_figure_triplet(figureHandle, folders.montage, names(typeIndex), ...
    parameters.export_dpi, "vector");
  clear cleanupFigure
end
end

function value = distribution_montage_title(plotType)
switch plotType
  case "domain"
    value = "2 deg orientation-domain ECD | number frequency; minPixel = 5";
  case "boundary"
    value = "Ti-Hex boundary misorientation | boundary-length weighted";
  otherwise
    value = "KAM distribution | first neighbour; 5 deg cutoff";
end
end

function draw_distribution_panel(axesHandle, state, plotType, edges, ...
    color, yMaximum, parameters, showLabels, showLegend)
centers = edges(1:end-1) + diff(edges) / 2;
switch plotType
  case "domain"
    frequency = state.domain_hist_percent;
  case "boundary"
    frequency = state.boundary_hist_percent;
  otherwise
    frequency = state.kam_hist_percent;
end
barHandle = bar(axesHandle, centers, frequency, 0.84, ...
  "FaceColor", color, "FaceAlpha", 0.74, ...
  "EdgeColor", [0.16 0.16 0.16], "LineWidth", 0.35);
hold(axesHandle, "on");
switch plotType
  case "domain"
    fineX = linspace(max(0.05, edges(1) + 0.05), edges(end), 600);
    fineY = 100 * median(diff(edges)) * lognormal_pdf(fineX, ...
      state.domain_fit_mu, state.domain_fit_sigma);
    fitHandle = plot(axesHandle, fineX, fineY, ...
      "Color", [0.08 0.08 0.08], "LineWidth", 1.35);
    if showLegend
      meanTextY = 0.70;
    else
      meanTextY = 0.92;
    end
    text(axesHandle, 0.97, meanTextY, sprintf("Mean ECD = %.2f um", ...
      state.domain_mean_um), "Units", "normalized", ...
      "HorizontalAlignment", "right", "VerticalAlignment", "top", ...
      "FontName", "Arial", "FontSize", 7.5, "FontWeight", "bold", ...
      "Color", color);
    if showLabels
      xlabel(axesHandle, "2 deg orientation-domain ECD (um)");
      ylabel(axesHandle, "Number frequency (%)");
    end
    if showLegend
      legend(axesHandle, [barHandle fitHandle], ...
        ["Histogram", "Descriptive lognormal fit"], ...
        "Location", "northeast", "FontName", "Arial", ...
        "FontSize", 6.5, "Box", "off");
    end
  case "boundary"
    xline(axesHandle, parameters.hagb_threshold_deg, "--", ...
      "15 deg", "Color", [0.15 0.15 0.15], "LineWidth", 0.9, ...
      "FontName", "Arial", "FontSize", 6.5, ...
      "LabelVerticalAlignment", "middle");
    if showLabels
      xlabel(axesHandle, "Boundary misorientation (deg)");
      ylabel(axesHandle, "Length frequency (%)");
    end
    text(axesHandle, 0.97, 0.92, sprintf("LAGB fraction = %.1f%%", ...
      100 * state.lagb_length_fraction), "Units", "normalized", ...
      "HorizontalAlignment", "right", "VerticalAlignment", "top", ...
      "FontName", "Arial", "FontSize", 7.5, "FontWeight", "bold", ...
      "Color", color);
  otherwise
    if showLabels
      xlabel(axesHandle, "KAM (deg)");
      ylabel(axesHandle, "Pixel frequency (%)");
    end
    text(axesHandle, 0.97, 0.92, sprintf("Mean KAM = %.3f deg", ...
      state.kam_mean_deg), "Units", "normalized", ...
      "HorizontalAlignment", "right", "VerticalAlignment", "top", ...
      "FontName", "Arial", "FontSize", 7.5, "FontWeight", "bold", ...
      "Color", color);
end
xlim(axesHandle, [edges(1) edges(end)]);
ylim(axesHandle, [0 yMaximum]);
set(axesHandle, "FontName", "Arial", "FontSize", 7, ...
  "TickDir", "out", "Box", "on", "LineWidth", 0.6, "Layer", "top");
end

function export_figure_triplet(figureHandle, folder, baseName, dpi, contentType)
drawnow;
exportgraphics(figureHandle, fullfile(folder, baseName + ".png"), ...
  "Resolution", dpi, "BackgroundColor", "white");
exportgraphics(figureHandle, fullfile(folder, baseName + ".tif"), ...
  "Resolution", dpi, "BackgroundColor", "white");
exportgraphics(figureHandle, fullfile(folder, baseName + ".pdf"), ...
  "ContentType", contentType, "Resolution", dpi, ...
  "BackgroundColor", "white");
end

function draw_ipf_inset(figureHandle, crystalSymmetry, position)
axesHandle = axes(figureHandle, "Units", "normalized", "Position", position);
colorKey = ipfHSVKey(crystalSymmetry);
colorKey.inversePoleFigureDirection = xvector;
plot(colorKey, "parent", axesHandle, "noTitle");
title(axesHandle, "IPF || AD", "FontName", "Arial", "FontSize", 6.5);
end

function draw_gb_inset_legend(axesHandle)
xLimits = xlim(axesHandle);
yLimits = ylim(axesHandle);
xRange = diff(xLimits);
yRange = diff(yLimits);
xLine = xLimits(1) + [0.64 0.73] * xRange;
xText = xLimits(1) + 0.75 * xRange;
yA = yLimits(1) + 0.94 * yRange;
yB = yLimits(1) + 0.88 * yRange;
line(axesHandle, xLine, [yA yA], ...
  "Color", [0.00 0.70 0.28], "LineWidth", 1.5);
text(axesHandle, xText, yA, "2-<15 deg", ...
  "VerticalAlignment", "middle", "FontName", "Arial", ...
  "FontSize", 6.5, "BackgroundColor", "white", "Margin", 0.5);
line(axesHandle, xLine, [yB yB], ...
  "Color", [0.03 0.03 0.03], "LineWidth", 1.5);
text(axesHandle, xText, yB, ">=15 deg", ...
  "VerticalAlignment", "middle", "FontName", "Arial", ...
  "FontSize", 6.5, "BackgroundColor", "white", "Margin", 0.5);
end

function draw_montage_gb_legend(figureHandle, position)
axesHandle = axes(figureHandle, "Units", "normalized", "Position", position);
axis(axesHandle, [0 1 0 1]);
axis(axesHandle, "off");
line(axesHandle, [0.03 0.18], [0.60 0.60], ...
  "Color", [0.00 0.70 0.28], "LineWidth", 1.8);
text(axesHandle, 0.21, 0.60, "LAGB: 2-<15 deg", ...
  "VerticalAlignment", "middle", "FontName", "Arial", "FontSize", 7);
line(axesHandle, [0.56 0.71], [0.60 0.60], ...
  "Color", [0.03 0.03 0.03], "LineWidth", 1.8);
text(axesHandle, 0.74, 0.60, "HAGB: >=15 deg", ...
  "VerticalAlignment", "middle", "FontName", "Arial", "FontSize", 7);
end

function draw_horizontal_color_scale(figureHandle, position, limits)
axesHandle = axes(figureHandle, "Units", "normalized", "Position", position);
imagesc(axesHandle, linspace(limits(1), limits(2), 256), 1, 1:256);
set(axesHandle, "YTick", [], "XTick", limits, "XLim", limits, ...
  "FontName", "Arial", "FontSize", 6.5, "Box", "on");
colormap(axesHandle, turbo(256));
xlabel(axesHandle, "KAM (deg)", "FontName", "Arial", "FontSize", 6.5);
end

function tag = state_tag(catalogRow)
diameter = replace(sprintf("%.2f", catalogRow.diameter_mm), ".", "p");
reduction = replace(sprintf("%.2f", ...
  catalogRow.cold_reduction_percent), ".", "p");
tag = diameter + "mm_" + reduction + "pct";
end

function limits = robust_limits(values)
values = values(isfinite(values));
limits = prctile(values, [1 99]);
if limits(1) >= limits(2)
  limits = [min(values) max(values)];
end
end

function edges = common_edges(values, width, lower, minimumUpper)
if nargin < 4
  minimumUpper = lower + width;
end
values = values(isfinite(values));
upper = width * ceil(max(values) / width);
upper = max(upper, minimumUpper);
edges = lower:width:upper;
if edges(end) < max(values)
  edges(end + 1) = edges(end) + width;
end
end

function value = common_ymax(frequency)
value = max(frequency, [], "all");
value = max(5, 5 * ceil(1.20 * value / 5));
end

function value = domain_ymax(state, edges)
fineX = linspace(max(0.05, edges(1) + 0.05), edges(end), 600);
peak = 0;
for stateIndex = 1:numel(state)
  fitted = 100 * median(diff(edges)) * lognormal_pdf(fineX, ...
    state(stateIndex).domain_fit_mu, state(stateIndex).domain_fit_sigma);
  peak = max([peak, state(stateIndex).domain_hist_percent, fitted]);
end
value = max(5, 5 * ceil(1.12 * peak / 5));
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

function format_map_axes(axesHandle)
axis(axesHandle, "image");
set(axesHandle, "YDir", "normal", "XTick", [], "YTick", [], ...
  "Box", "on", "LineWidth", 0.45, "FontName", "Arial");
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
  "Color", "black", "LineWidth", 1.6);
text(axesHandle, xStart + lengthUm / 2, yPosition + 0.018 * yRange, ...
  sprintf("%d um", lengthUm), "HorizontalAlignment", "center", ...
  "VerticalAlignment", "bottom", "FontName", "Arial", ...
  "FontSize", 6, "BackgroundColor", "white", "Margin", 0.5);
end

function percentages = weighted_histogram(values, weights, edges)
values = double(values(:));
weights = double(weights(:));
valid = isfinite(values) & isfinite(weights) & weights > 0;
values = values(valid);
weights = weights(valid);
bin = discretize(values, edges);
assert(all(isfinite(bin)), "Common bins do not span all values.");
percentages = 100 * accumarray(bin, weights, ...
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
values = double(values);
density = exp(-0.5 * ((log(values) - mu) / sigma).^2) ./ ...
  (values * sigma * sqrt(2 * pi));
end

function q = weighted_quantile(values, weights, probability)
values = double(values(:));
weights = double(weights(:));
valid = isfinite(values) & isfinite(weights) & weights > 0;
values = values(valid);
weights = weights(valid);
[values, order] = sort(values);
weights = weights(order);
index = find(cumsum(weights) >= probability * sum(weights), 1);
q = values(index);
end

function summary = build_summary(state, catalog, parameters)
sample = catalog.sample;
diameter_mm = catalog.diameter_mm;
cold_reduction_percent = catalog.cold_reduction_percent;
variant = catalog.variant;
input_path = catalog.input_path;
orientation_domain_count = reshape([state.domain_count], [], 1);
orientation_domain_ecd_number_mean_um = reshape( ...
  [state.domain_mean_um], [], 1);
orientation_domain_ecd_number_median_um = reshape( ...
  [state.domain_median_um], [], 1);
orientation_domain_ecd_area_weighted_mean_um = reshape( ...
  [state.domain_area_weighted_mean_um], [], 1);
orientation_domain_ecd_area_weighted_median_um = reshape( ...
  [state.domain_area_weighted_median_um], [], 1);
lagb_length_fraction = reshape([state.lagb_length_fraction], [], 1);
hagb_length_fraction = reshape([state.hagb_length_fraction], [], 1);
kam_mean_deg = reshape([state.kam_mean_deg], [], 1);
kam_median_deg = reshape([state.kam_median_deg], [], 1);
kam_p90_deg = reshape([state.kam_p90_deg], [], 1);
boundary_detection_deg = repmat(parameters.boundary_detection_deg, 6, 1);
hagb_threshold_deg = repmat(parameters.hagb_threshold_deg, 6, 1);
min_domain_pixels = repmat(parameters.min_domain_pixels, 6, 1);
kam_order = repmat(parameters.kam_order, 6, 1);
kam_cutoff_deg = repmat(parameters.kam_cutoff_deg, 6, 1);
summary = table(sample, diameter_mm, cold_reduction_percent, variant, ...
  input_path, orientation_domain_count, ...
  orientation_domain_ecd_number_mean_um, ...
  orientation_domain_ecd_number_median_um, ...
  orientation_domain_ecd_area_weighted_mean_um, ...
  orientation_domain_ecd_area_weighted_median_um, ...
  lagb_length_fraction, hagb_length_fraction, kam_mean_deg, ...
  kam_median_deg, kam_p90_deg, boundary_detection_deg, ...
  hagb_threshold_deg, min_domain_pixels, kam_order, kam_cutoff_deg);
end

function output = build_distribution_table(state, catalog, domainEdges, ...
    boundaryEdges, kamEdges)
plot_type = strings(0, 1);
sample = strings(0, 1);
cold_reduction_percent = zeros(0, 1);
bin_lower = zeros(0, 1);
bin_upper = zeros(0, 1);
bin_center = zeros(0, 1);
frequency_percent = zeros(0, 1);
for stateIndex = 1:6
  [plot_type, sample, cold_reduction_percent, bin_lower, bin_upper, ...
    bin_center, frequency_percent] = append_distribution(plot_type, ...
    sample, cold_reduction_percent, bin_lower, bin_upper, bin_center, ...
    frequency_percent, "orientation_domain_ecd_number", ...
    catalog(stateIndex, :), domainEdges, ...
    state(stateIndex).domain_hist_percent);
  [plot_type, sample, cold_reduction_percent, bin_lower, bin_upper, ...
    bin_center, frequency_percent] = append_distribution(plot_type, ...
    sample, cold_reduction_percent, bin_lower, bin_upper, bin_center, ...
    frequency_percent, "boundary_misorientation_length_weighted", ...
    catalog(stateIndex, :), boundaryEdges, ...
    state(stateIndex).boundary_hist_percent);
  [plot_type, sample, cold_reduction_percent, bin_lower, bin_upper, ...
    bin_center, frequency_percent] = append_distribution(plot_type, ...
    sample, cold_reduction_percent, bin_lower, bin_upper, bin_center, ...
    frequency_percent, "kam_pixel_frequency", catalog(stateIndex, :), ...
    kamEdges, state(stateIndex).kam_hist_percent);
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

function validate_outputs(state, domainEdges, boundaryEdges, kamEdges, folders)
for stateIndex = 1:6
  assert(abs(sum(state(stateIndex).domain_hist_percent) - 100) < 1e-8);
  assert(abs(sum(state(stateIndex).boundary_hist_percent) - 100) < 1e-8);
  assert(abs(sum(state(stateIndex).kam_hist_percent) - 100) < 1e-8);
end
assert(domainEdges(1) == 0 && boundaryEdges(1) == 0 && kamEdges(1) == 0);
assert(numel(dir(fullfile(folders.ipf, "*.png"))) == 6);
assert(numel(dir(fullfile(folders.gb, "*.png"))) == 6);
assert(numel(dir(fullfile(folders.kam, "*.png"))) == 6);
assert(numel(dir(fullfile(folders.domain, "*.png"))) == 6);
assert(numel(dir(fullfile(folders.boundary, "*.png"))) == 6);
assert(numel(dir(fullfile(folders.kamdist, "*.png"))) == 6);
assert(numel(dir(fullfile(folders.montage, "*.png"))) == 6);
fprintf("COMPONENT_GALLERY validation passed: 36 panels + 6 montages.\n");
end
