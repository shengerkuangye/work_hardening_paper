function summary = generate_fig7_texture_gallery(scanRoot, outputDir)
%GENERATE_FIG7_TEXTURE_GALLERY Six-state PF/IPF gallery after reference Fig. 7.
% Each state contains {0001}, {10-10}, and {11-20} pole figures plus
% inverse pole figures for AD(X), TD/RD(Y), and ND(Z). All pole figures
% share one MRD scale, and all inverse pole figures share another MRD
% scale, so the six deformation states can be compared directly.

arguments
  scanRoot (1,1) string
  outputDir (1,1) string
end

assert(isfolder(scanRoot), "EBSD scan folder not found: %s", scanRoot);
assert(~isempty(which("EBSD")), "MTEX must be loaded before this function.");
if ~isfolder(outputDir), mkdir(outputDir); end
stateDir = fullfile(outputDir, "01_state_figures");
montageDir = fullfile(outputDir, "montages");
if ~isfolder(stateDir), mkdir(stateDir); end
if ~isfolder(montageDir), mkdir(montageDir); end

parameters = struct( ...
  "kernel_halfwidth_deg", 5, ...
  "plot_resolution_deg", 5, ...
  "export_dpi", 600, ...
  "coordinate_x", "AD", ...
  "coordinate_y", "TD/RD", ...
  "coordinate_z", "ND", ...
  "pole_figures", "{0001}; {10-10}; {11-20}", ...
  "inverse_pole_directions", "AD(X); TD/RD(Y); ND(Z)");

catalog = comprehensive_ebsd_catalog(scanRoot);
catalog = catalog(catalog.variant == "raw", :);
[~, order] = sort(catalog.cold_reduction_percent);
catalog = catalog(order, :);
assert(height(catalog) == 6, "Expected six registered raw EBSD states.");

kernel = SO3DeLaValleePoussinKernel( ...
  "halfwidth", parameters.kernel_halfwidth_deg * degree);
odfs = cell(6, 1);
orientationCount = zeros(6, 1);
pfAutoMax = zeros(6, 1);
ipfAutoMax = zeros(6, 1);
crystalSymmetry = [];

pfAnnotations = @(varargin) text( ...
  [vector3d.X, vector3d.Y, vector3d.Z], ...
  {"", "", ""}, "Visible", "off", ...
  "FontName", "Arial", "FontSize", 4.5, "FontWeight", "normal", ...
  "Margin", 0.25, "tag", "axesLabels", varargin{:});
setMTEXpref("pfAnnotations", pfAnnotations);

for stateIndex = 1:6
  fprintf("FIG7_TEXTURE ODF sample=%s reduction=%.2f%%\n", ...
    catalog.sample(stateIndex), catalog.cold_reduction_percent(stateIndex));
  [ebsdFull, ~] = load_comprehensive_ebsd_scan(catalog(stateIndex, :));
  tiEbsd = ebsdFull("Ti-Hex");
  assert(~isempty(tiEbsd), "Ti-Hex phase is absent.");
  orientationCount(stateIndex) = length(tiEbsd);
  if isempty(crystalSymmetry)
    crystalSymmetry = tiEbsd.CS;
  end
  rbfOdf = calcDensity(tiEbsd.orientations, "kernel", kernel, ...
    "weights", ones(length(tiEbsd), 1), "silent");
  rbfOdf = normalize_positive_mean_density(rbfOdf);
  odfs{stateIndex} = SO3FunHarmonic(rbfOdf, ...
    "bandwidth", kernel.bandwidth);
  odfs{stateIndex} = normalize_positive_mean_density(odfs{stateIndex});

  poles = texture_poles(tiEbsd.CS);
  specimenDirections = [xvector, yvector, zvector];
  pfAutoMax(stateIndex) = automatic_plot_maximum(odfs{stateIndex}, ...
    "pf", poles, specimenDirections, parameters.plot_resolution_deg);
  ipfAutoMax(stateIndex) = automatic_plot_maximum(odfs{stateIndex}, ...
    "ipf", poles, specimenDirections, parameters.plot_resolution_deg);
  clear ebsdFull tiEbsd rbfOdf poles specimenDirections
end

globalPfMaximum = ceil(max(pfAutoMax));
globalIpfMaximum = ceil(max(ipfAutoMax));
globalPfMaximum = max(globalPfMaximum, 1);
globalIpfMaximum = max(globalIpfMaximum, 1);
fprintf("FIG7_TEXTURE common PF max=%.2f MRD; IPF max=%.2f MRD\n", ...
  globalPfMaximum, globalIpfMaximum);

tempRoot = string(tempname);
mkdir(tempRoot);
cleanupTemp = onCleanup(@() remove_temp_directory(tempRoot));
pfRowPaths = strings(6, 1);
ipfRowPaths = strings(6, 1);
for stateIndex = 1:6
  pfRowPaths(stateIndex) = fullfile(tempRoot, ...
    sprintf("pf_row_%d.png", stateIndex));
  ipfRowPaths(stateIndex) = fullfile(tempRoot, ...
    sprintf("ipf_row_%d.png", stateIndex));
  poles = texture_poles(crystalSymmetry);
  specimenDirections = [xvector, yvector, zvector];
  render_texture_row(odfs{stateIndex}, "pf", poles, ...
    specimenDirections, globalPfMaximum, ...
    parameters.plot_resolution_deg, pfRowPaths(stateIndex));
  render_texture_row(odfs{stateIndex}, "ipf", poles, ...
    specimenDirections, globalIpfMaximum, ...
    parameters.plot_resolution_deg, ipfRowPaths(stateIndex));
end
crop_rows_to_common_margins(pfRowPaths, 30);
crop_rows_to_common_margins(ipfRowPaths, 30);

statePaths = strings(6, 1);
for stateIndex = 1:6
  tag = state_tag(catalog(stateIndex, :));
  statePaths(stateIndex) = fullfile(stateDir, tag + "_pf_ipf_texture.png");
  render_state_figure(pfRowPaths(stateIndex), ipfRowPaths(stateIndex), ...
    catalog(stateIndex, :), globalPfMaximum, globalIpfMaximum, ...
    parameters.export_dpi, statePaths(stateIndex));
end

montageBase = fullfile(montageDir, "fig7_texture_six_state_montage");
render_montage(statePaths, montageBase, parameters.export_dpi);

sample = catalog.sample;
diameter_mm = catalog.diameter_mm;
cold_reduction_percent = catalog.cold_reduction_percent;
variant = catalog.variant;
input_path = catalog.input_path;
valid_ti_hex_orientation_count = orientationCount;
kernel_halfwidth_deg = repmat(parameters.kernel_halfwidth_deg, 6, 1);
plot_resolution_deg = repmat(parameters.plot_resolution_deg, 6, 1);
pf_auto_maximum_mrd = pfAutoMax;
ipf_auto_maximum_mrd = ipfAutoMax;
global_pf_color_maximum_mrd = repmat(globalPfMaximum, 6, 1);
global_ipf_color_maximum_mrd = repmat(globalIpfMaximum, 6, 1);
coordinate_x = repmat("AD", 6, 1);
coordinate_y = repmat("TD/RD", 6, 1);
coordinate_z = repmat("ND", 6, 1);
summary = table(sample, diameter_mm, cold_reduction_percent, variant, ...
  input_path, valid_ti_hex_orientation_count, kernel_halfwidth_deg, ...
  plot_resolution_deg, pf_auto_maximum_mrd, ipf_auto_maximum_mrd, ...
  global_pf_color_maximum_mrd, global_ipf_color_maximum_mrd, ...
  coordinate_x, coordinate_y, coordinate_z);
summary.output_file = statePaths;
writetable(summary, fullfile(outputDir, "fig7_texture_summary.csv"));
writetable(struct2table(parameters), ...
  fullfile(outputDir, "fig7_texture_parameters.csv"));

assert(numel(dir(fullfile(stateDir, "*_pf_ipf_texture.png"))) == 6);
assert(isfile(montageBase + ".png"));
assert(isfile(montageBase + ".tif"));
assert(isfile(montageBase + ".pdf"));
fprintf("FIG7_TEXTURE validation passed: 6 state figures + 1 montage.\n");
clear cleanupTemp
end

function poles = texture_poles(crystalSymmetry)
poles = [ ...
  Miller(0, 0, 0, 1, crystalSymmetry), ...
  Miller(1, 0, -1, 0, crystalSymmetry), ...
  Miller(1, 1, -2, 0, crystalSymmetry)];
end

function maximum = automatic_plot_maximum(odf, plotType, poles, ...
    specimenDirections, resolutionDeg)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Units", "pixels", "Position", [50 50 1100 390]);
cleanupFigure = onCleanup(@() close(figureHandle));
if plotType == "pf"
  plotPDF(odf, poles, "antipodal", "contourf", "silent", ...
    "resolution", resolutionDeg * degree, "layout", [1 3]);
else
  plotIPDF(odf, specimenDirections, "antipodal", "contourf", ...
    "silent", "resolution", resolutionDeg * degree, "layout", [1 3]);
end
drawnow;
axesHandles = spherical_axes(figureHandle);
assert(numel(axesHandles) == 3, "MTEX did not render three texture axes.");
limits = vertcat(axesHandles.CLim);
maximum = max(limits(:, 2));
assert(isfinite(maximum) && maximum > 0);
clear cleanupFigure
end

function render_texture_row(odf, plotType, poles, specimenDirections, ...
    colorMaximum, resolutionDeg, outputPath)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Units", "inches", "Position", [0.25 0.25 8.2 2.65]);
cleanupFigure = onCleanup(@() close(figureHandle));
if plotType == "pf"
  plotPDF(odf, poles, "antipodal", "contourf", "silent", ...
    "resolution", resolutionDeg * degree, "layout", [1 3]);
  panelTitles = ["{0001}", "{10-10}", "{11-20}"];
else
  plotIPDF(odf, specimenDirections, "antipodal", "contourf", ...
    "silent", "resolution", resolutionDeg * degree, "layout", [1 3]);
  panelTitles = ["AD (X)", "TD/RD (Y)", "ND (Z)"];
end
setColorRange([0 colorMaximum], "current");
mtexColorMap turbo;
axesHandles = sort_axes_left_to_right(spherical_axes(figureHandle));
for panelIndex = 1:3
  titleHandle = title(axesHandles(panelIndex), panelTitles(panelIndex), ...
    "FontName", "Arial", "FontSize", 8, "FontWeight", "bold", ...
    "Interpreter", "none");
  if plotType == "pf"
    % MTEX places the specimen-direction labels at the pole-figure rim.
    % Lift the Miller-index title above them and keep the coordinate labels
    % compact so neither item obscures the density contours.
    titleHandle.Units = "normalized";
    titlePosition = titleHandle.Position;
    titlePosition(2) = 1.08;
    titleHandle.Position = titlePosition;
    delete(findall(figureHandle, "Tag", "axesLabels"));
  end
end
set(figureHandle, "Units", "inches", ...
  "Position", [0.25 0.25 8.2 2.65]);
drawnow;
exportgraphics(figureHandle, outputPath, "Resolution", 600, ...
  "BackgroundColor", "white");
clear cleanupFigure
end

function render_state_figure(pfPath, ipfPath, catalogRow, ...
    pfMaximum, ipfMaximum, dpi, outputPath)
pfImage = imread(pfPath);
ipfImage = imread(ipfPath);
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Units", "inches", "Position", [0.25 0.25 8.6 5.9]);
cleanupFigure = onCleanup(@() close(figureHandle));

annotation(figureHandle, "textbox", [0.02 0.952 0.96 0.038], ...
  "String", sprintf("%.2f%% cold reduction | %.2f mm", ...
  catalogRow.cold_reduction_percent, catalogRow.diameter_mm), ...
  "LineStyle", "none", "HorizontalAlignment", "center", ...
  "VerticalAlignment", "middle", "FontName", "Arial", ...
  "FontSize", 10, "FontWeight", "bold");

pfAxes = axes(figureHandle, "Units", "normalized", ...
  "Position", [0.040 0.525 0.85 0.395]);
image(pfAxes, pfImage);
axis(pfAxes, "image"); axis(pfAxes, "off");
ipfAxes = axes(figureHandle, "Units", "normalized", ...
  "Position", [0.040 0.070 0.85 0.395]);
image(ipfAxes, ipfImage);
axis(ipfAxes, "image"); axis(ipfAxes, "off");

annotation(figureHandle, "textbox", [0.008 0.850 0.10 0.035], ...
  "String", "PF", "LineStyle", "none", "FontName", "Arial", ...
  "FontSize", 8, "FontWeight", "bold");
annotation(figureHandle, "textbox", [0.008 0.390 0.10 0.035], ...
  "String", "IPF", "LineStyle", "none", "FontName", "Arial", ...
  "FontSize", 8, "FontWeight", "bold");
draw_vertical_color_scale(figureHandle, [0.915 0.565 0.018 0.295], ...
  pfMaximum);
draw_vertical_color_scale(figureHandle, [0.915 0.110 0.018 0.295], ...
  ipfMaximum);
annotation(figureHandle, "textbox", [0.02 0.008 0.96 0.030], ...
  "String", sprintf("Raw Ti-Hex EBSD | PF coordinates: X=AD, Y=TD/RD, center=ND | common scales: PF 0-%.0f MRD; IPF 0-%.0f MRD", ...
  pfMaximum, ipfMaximum), ...
  "LineStyle", "none", "HorizontalAlignment", "center", ...
  "FontName", "Arial", "FontSize", 7);
drawnow;
exportgraphics(figureHandle, outputPath, "Resolution", dpi, ...
  "BackgroundColor", "white");
clear cleanupFigure
end

function draw_vertical_color_scale(figureHandle, position, maximum)
axesHandle = axes(figureHandle, "Units", "normalized", ...
  "Position", position);
imagesc(axesHandle, 1, linspace(0, maximum, 256), (1:256).');
set(axesHandle, "YDir", "normal", "XTick", [], ...
  "YTick", [0 maximum], "YLim", [0 maximum], ...
  "YAxisLocation", "right", "FontName", "Arial", ...
  "FontSize", 6.5, "Box", "on", "LineWidth", 0.5);
colormap(axesHandle, turbo(256));
ylabel(axesHandle, "m.r.d.", "FontName", "Arial", "FontSize", 6.5);
end

function render_montage(statePaths, montageBase, dpi)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Units", "inches", "Position", [0.25 0.25 15.0 7.8]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 2, 3, "Padding", "compact", ...
  "TileSpacing", "compact");
for stateIndex = 1:6
  axesHandle = nexttile(layout);
  stateImage = imread(statePaths(stateIndex));
  image(axesHandle, stateImage);
  axis(axesHandle, "image");
  axis(axesHandle, "off");
  text(axesHandle, 0.012, 0.985, sprintf("(%c)", ...
    char('a' + stateIndex - 1)), "Units", "normalized", ...
    "HorizontalAlignment", "left", "VerticalAlignment", "top", ...
    "FontName", "Arial", "FontSize", 8, "FontWeight", "bold", ...
    "BackgroundColor", "white", "Margin", 1);
end
title(layout, ["Texture analysis after reference Fig. 7 | " ...
  "pole figures and inverse pole figures"], ...
  "FontName", "Arial", "FontSize", 11, "FontWeight", "bold");
drawnow;
exportgraphics(figureHandle, montageBase + ".png", ...
  "Resolution", dpi, "BackgroundColor", "white");
exportgraphics(figureHandle, montageBase + ".tif", ...
  "Resolution", dpi, "BackgroundColor", "white");
exportgraphics(figureHandle, montageBase + ".pdf", ...
  "ContentType", "image", "Resolution", dpi, ...
  "BackgroundColor", "white");
clear cleanupFigure
end

function handles = spherical_axes(figureHandle)
handles = findall(figureHandle, "Type", "axes");
handles = handles(arrayfun(@(axisHandle) ...
  isappdata(axisHandle, "sphericalPlot"), handles));
end

function handles = sort_axes_left_to_right(handles)
positions = vertcat(handles.Position);
[~, order] = sort(positions(:, 1));
handles = handles(order);
end

function crop_rows_to_common_margins(paths, padding)
images = cellfun(@imread, cellstr(paths), "UniformOutput", false);
sizes = cellfun(@size, images, "UniformOutput", false);
assert(all(cellfun(@(imageSize) isequal(imageSize, sizes{1}), sizes)));
referenceSize = sizes{1};
unionNonwhite = false(referenceSize(1), referenceSize(2));
for imageIndex = 1:numel(images)
  imageData = images{imageIndex};
  if ismatrix(imageData)
    nonwhite = imageData < 250;
  else
    nonwhite = any(imageData < 250, 3);
  end
  unionNonwhite = unionNonwhite | nonwhite;
end
[contentRows, contentColumns] = find(unionNonwhite);
assert(~isempty(contentRows), "Texture rows are blank.");
firstRow = max(1, min(contentRows) - padding);
lastRow = min(referenceSize(1), max(contentRows) + padding);
firstColumn = max(1, min(contentColumns) - padding);
lastColumn = min(referenceSize(2), max(contentColumns) + padding);
for imageIndex = 1:numel(images)
  imageData = images{imageIndex};
  imwrite(imageData(firstRow:lastRow, firstColumn:lastColumn, :), ...
    paths(imageIndex));
end
end

function tag = state_tag(catalogRow)
diameter = replace(sprintf("%.2f", catalogRow.diameter_mm), ".", "p");
reduction = replace(sprintf("%.2f", ...
  catalogRow.cold_reduction_percent), ".", "p");
tag = diameter + "mm_" + reduction + "pct";
end

function remove_temp_directory(path)
if isfolder(path)
  rmdir(path, "s");
end
end
