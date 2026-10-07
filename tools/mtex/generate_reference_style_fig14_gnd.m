function outputs = generate_reference_style_fig14_gnd(scanRoot, outputRoot)
%GENERATE_REFERENCE_STYLE_FIG14_GND Create six-state KAM-derived GND figure.
% The layout adapts Xia et al. (2026) Fig. 14 to all six project states.
% KAM is converted pixelwise to an estimated GND density. This scalar
% estimate is an orientation-gradient proxy and is not total dislocation
% density.

arguments
  scanRoot (1,1) string
  outputRoot (1,1) string
end

assert(isfolder(scanRoot), "EBSD scan folder not found: %s", scanRoot);
assert(~isempty(which("EBSD")), ...
  "MTEX is not loaded in the current MATLAB session.");
if ~isfolder(outputRoot)
  mkdir(outputRoot);
end

parameters = struct( ...
  "input_variant", "raw", ...
  "grain_detection_deg", 2, ...
  "kam_order", 1, ...
  "kam_cutoff_deg", 5, ...
  "burgers_vector_nm", 0.295, ...
  "gnd_map_max_1e14_m2", 8, ...
  "gnd_hist_max_1e14_m2", 12, ...
  "gnd_hist_bin_width_1e14_m2", 0.1, ...
  "scale_bar_um", 100, ...
  "export_dpi", 600);

catalog = comprehensive_ebsd_catalog(scanRoot);
catalog = catalog(catalog.variant == parameters.input_variant, :);
[~, order] = sort(catalog.cold_reduction_percent);
catalog = catalog(order, :);
assert(height(catalog) == 6, "Expected six registered raw EBSD states.");

state = repmat(empty_state(), height(catalog), 1);
for stateIndex = 1:height(catalog)
  fprintf("REFERENCE_FIG14 sample=%s reduction=%.2f%%\n", ...
    catalog.sample(stateIndex), catalog.cold_reduction_percent(stateIndex));
  state(stateIndex) = calculate_state(catalog(stateIndex, :), parameters);
end

histEdges = 0:parameters.gnd_hist_bin_width_1e14_m2: ...
  parameters.gnd_hist_max_1e14_m2;
for stateIndex = 1:numel(state)
  state(stateIndex).hist_percent = pixel_histogram( ...
    state(stateIndex).gnd_values_1e14_m2, histEdges);
end

baseName = "reference_style_fig14_gnd_six_diameters";
pngPath = fullfile(outputRoot, baseName + ".png");
tifPath = fullfile(outputRoot, baseName + ".tif");
pdfPath = fullfile(outputRoot, baseName + ".pdf");
render_figure(state, catalog, histEdges, parameters, ...
  pngPath, tifPath, pdfPath);

summary = build_summary(state, catalog, parameters);
summaryPath = fullfile(outputRoot, ...
  "reference_style_fig14_gnd_summary.csv");
distributionPath = fullfile(outputRoot, ...
  "reference_style_fig14_gnd_distribution.csv");
parameterPath = fullfile(outputRoot, ...
  "reference_style_fig14_gnd_parameters.csv");
writetable(summary, summaryPath);
writetable(build_distribution_table(state, catalog, histEdges), ...
  distributionPath);
writetable(build_parameter_table(parameters, state), parameterPath);
write_readme(outputRoot, parameters);

validate_outputs(state, summary, histEdges, ...
  [pngPath; tifPath; pdfPath; summaryPath; distributionPath; parameterPath]);
outputs = struct("png", pngPath, "tif", tifPath, "pdf", pdfPath, ...
  "summary", summaryPath, "distribution", distributionPath, ...
  "parameters", parameterPath);
end

function result = empty_state()
result = struct( ...
  "x", [], "y", [], "gnd_image_1e14_m2", [], ...
  "gnd_values_1e14_m2", [], "hist_percent", [], ...
  "step_um", NaN, ...
  "kam_mean_deg", NaN, "kam_median_deg", NaN, "kam_p90_deg", NaN, ...
  "gnd_mean_1e14_m2", NaN, "gnd_median_1e14_m2", NaN, ...
  "gnd_p90_1e14_m2", NaN, "map_saturation_percent", NaN, ...
  "valid_pixel_count", 0);
end

function result = calculate_state(catalogRow, parameters)
[ebsdFull, ~] = load_comprehensive_ebsd_scan(catalogRow);
assert(~isempty(ebsdFull("Ti-Hex")), "Ti-Hex phase is absent.");

[~, kamGrainId] = calcGrains(ebsdFull, "unitCell", ...
  "threshold", parameters.grain_detection_deg * degree);
ebsdFull.grainId = kamGrainId;
ebsdGrid = ebsdFull.gridify;
kam = ebsdGrid.KAM("order", parameters.kam_order, ...
  "threshold", parameters.kam_cutoff_deg * degree);
kamDeg = reshape(double(kam / degree), [], 1);

tiGrid = ebsdGrid("Ti-Hex");
phaseId = unique(double(tiGrid.phaseId));
assert(isscalar(phaseId), "Expected one Ti-Hex phase identifier.");
tiMask = reshape(double(ebsdGrid.phaseId), [], 1) == phaseId;
kamDeg(~tiMask) = NaN;

[xValues, yValues, nativeIndices] = native_grid_indices(ebsdGrid);
stepX = median(diff(xValues));
stepY = median(diff(yValues));
assert(isfinite(stepX) && isfinite(stepY) && stepX > 0 && stepY > 0);
assert(abs(stepX - stepY) <= 1e-8 * max(stepX, stepY), ...
  "GND conversion requires a square EBSD grid.");
stepUm = mean([stepX stepY]);

bMeters = parameters.burgers_vector_nm * 1e-9;
stepMeters = stepUm * 1e-6;
gndM2 = 2 * deg2rad(kamDeg) / (stepMeters * bMeters);
gndScaled = gndM2 / 1e14;

valid = isfinite(gndScaled);
values = gndScaled(valid);
assert(~isempty(values), "No finite Ti-Hex KAM pixels were found.");

result = empty_state();
result.x = xValues;
result.y = yValues;
result.gnd_image_1e14_m2 = scalar_image(gndScaled, nativeIndices, ...
  numel(yValues), numel(xValues));
result.gnd_values_1e14_m2 = values;
result.step_um = stepUm;
result.kam_mean_deg = mean(kamDeg(valid));
result.kam_median_deg = median(kamDeg(valid));
result.kam_p90_deg = prctile(kamDeg(valid), 90);
result.gnd_mean_1e14_m2 = mean(values);
result.gnd_median_1e14_m2 = median(values);
result.gnd_p90_1e14_m2 = prctile(values, 90);
result.map_saturation_percent = 100 * mean( ...
  values > parameters.gnd_map_max_1e14_m2);
result.valid_pixel_count = numel(values);
end

function render_figure(state, catalog, histEdges, parameters, ...
    pngPath, tifPath, pdfPath)
colors = [ ...
  0.12 0.35 0.80
  0.00 0.58 0.74
  0.12 0.66 0.42
  0.84 0.64 0.08
  0.92 0.38 0.10
  0.76 0.12 0.18];

figureHandle = figure("Visible", "off", "Color", "white", ...
  "Units", "inches", "Position", [0.2 0.2 12.4 10.2]);
cleanupFigure = onCleanup(@() close(figureHandle));

mapLeft = [0.035 0.305];
mapWidth = 0.235;
mapHeight = 0.235;
histLeft = 0.585;
histWidth = 0.385;
histHeight = 0.105;
rowBottom = [0.675 0.375 0.075];

for pairIndex = 1:3
  firstState = 2 * pairIndex - 1;
  secondState = firstState + 1;
  panelLetters = char('a' + (pairIndex - 1) * 3 + (0:2));

  for mapIndex = 1:2
    stateIndex = firstState + mapIndex - 1;
    axesHandle = axes(figureHandle, "Units", "normalized", ...
      "Position", [mapLeft(mapIndex) rowBottom(pairIndex) ...
      mapWidth mapHeight]);
    draw_gnd_map(axesHandle, state(stateIndex), parameters);
    title(axesHandle, sprintf("(%c) %s: %.2f mm, %.2f%% reduction", ...
      panelLetters(mapIndex), catalog.sample(stateIndex), ...
      catalog.diameter_mm(stateIndex), ...
      catalog.cold_reduction_percent(stateIndex)), ...
      "FontName", "Arial", "FontSize", 8.2, "FontWeight", "bold");
  end

  pairFrequency = [state(firstState).hist_percent; ...
    state(secondState).hist_percent];
  pairYMax = max(5, 5 * ceil(1.15 * max(pairFrequency, [], "all") / 5));
  for histIndex = 1:2
    stateIndex = firstState + histIndex - 1;
    histBottom = rowBottom(pairIndex) + ...
      (2 - histIndex) * (histHeight + 0.014);
    axesHandle = axes(figureHandle, "Units", "normalized", ...
      "Position", [histLeft histBottom histWidth histHeight]);
    draw_histogram(axesHandle, state(stateIndex), catalog(stateIndex, :), ...
      histEdges, colors(stateIndex, :), pairYMax, histIndex == 2);
    if histIndex == 1
      title(axesHandle, sprintf("(%c) GND distributions", panelLetters(3)), ...
        "FontName", "Arial", "FontSize", 8.2, "FontWeight", "bold");
    end
  end
end

annotation(figureHandle, "textbox", [0.01 0.965 0.98 0.025], ...
  "String", "KAM-derived GND density maps and distributions", ...
  "LineStyle", "none", "HorizontalAlignment", "center", ...
  "FontName", "Arial", "FontSize", 11, "FontWeight", "bold");
annotation(figureHandle, "textbox", [0.01 0.94 0.98 0.022], ...
  "String", sprintf("Raw EBSD; first-neighbour KAM; %.0f deg cutoff; " + ...
  "step = %.1f um; b = %.3f nm; common map scale", ...
  parameters.kam_cutoff_deg, state(1).step_um, ...
  parameters.burgers_vector_nm), "LineStyle", "none", ...
  "HorizontalAlignment", "center", "FontName", "Arial", ...
  "FontSize", 7.5);
draw_color_scale(figureHandle, [0.25 0.036 0.50 0.018], parameters);

drawnow;
exportgraphics(figureHandle, pngPath, "Resolution", parameters.export_dpi, ...
  "BackgroundColor", "white");
exportgraphics(figureHandle, tifPath, "Resolution", parameters.export_dpi, ...
  "BackgroundColor", "white");
exportgraphics(figureHandle, pdfPath, "ContentType", "image", ...
  "Resolution", parameters.export_dpi, "BackgroundColor", "white");
clear cleanupFigure
end

function draw_gnd_map(axesHandle, current, parameters)
imageHandle = imagesc(axesHandle, current.x, current.y, ...
  current.gnd_image_1e14_m2);
set(imageHandle, "AlphaData", isfinite(current.gnd_image_1e14_m2));
axis(axesHandle, "image");
set(axesHandle, "YDir", "normal", "XTick", [], "YTick", [], ...
  "Box", "on", "LineWidth", 0.5, "Color", [0.78 0.78 0.78]);
hold(axesHandle, "on");
colormap(axesHandle, turbo(256));
clim(axesHandle, [0 parameters.gnd_map_max_1e14_m2]);
draw_scale_bar(axesHandle, current.x, current.y, parameters.scale_bar_um);
text(axesHandle, 0.98, 0.965, sprintf("mean = %.2f x 10^{14} m^{-2}\n" + ...
  "KAM = %.3f deg", current.gnd_mean_1e14_m2, current.kam_mean_deg), ...
  "Units", "normalized", "HorizontalAlignment", "right", ...
  "VerticalAlignment", "top", "FontName", "Arial", "FontSize", 6.4, ...
  "BackgroundColor", "white", "Margin", 1);
end

function draw_histogram(axesHandle, current, catalogRow, edges, color, ...
    yMaximum, showXLabel)
centers = edges(1:end-1) + diff(edges) / 2;
bar(axesHandle, centers, current.hist_percent, 1, ...
  "FaceColor", color, "FaceAlpha", 0.78, ...
  "EdgeColor", color, "LineWidth", 0.15);
xlim(axesHandle, [edges(1) edges(end)]);
ylim(axesHandle, [0 yMaximum]);
set(axesHandle, "FontName", "Arial", "FontSize", 6.5, ...
  "TickDir", "out", "Box", "on", "LineWidth", 0.55, "Layer", "top");
ylabel(axesHandle, "Frequency (%)", "FontName", "Arial", "FontSize", 6.5);
if showXLabel
  xlabel(axesHandle, "GND density (10^{14} m^{-2})", ...
    "FontName", "Arial", "FontSize", 6.5);
else
  set(axesHandle, "XTickLabel", []);
end
text(axesHandle, 0.025, 0.88, sprintf("%s | mean = %.2f x 10^{14} m^{-2}", ...
  catalogRow.sample, current.gnd_mean_1e14_m2), ...
  "Units", "normalized", "HorizontalAlignment", "left", ...
  "VerticalAlignment", "top", "FontName", "Arial", ...
  "FontSize", 6.7, "FontWeight", "bold", "Color", color);
end

function draw_color_scale(figureHandle, position, parameters)
axesHandle = axes(figureHandle, "Units", "normalized", "Position", position);
limits = [0 parameters.gnd_map_max_1e14_m2];
imagesc(axesHandle, linspace(limits(1), limits(2), 256), 1, 1:256);
set(axesHandle, "YTick", [], "XTick", 0:2:limits(2), ...
  "XLim", limits, "FontName", "Arial", "FontSize", 6.5, "Box", "on");
colormap(axesHandle, turbo(256));
xlabel(axesHandle, "Estimated GND density (10^{14} m^{-2})", ...
  "FontName", "Arial", "FontSize", 6.7);
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
  "FontSize", 5.8, "BackgroundColor", "white", "Margin", 0.5);
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

function percentages = pixel_histogram(values, edges)
values = double(values(:));
values = values(isfinite(values));
assert(~isempty(values));
assert(min(values) >= edges(1));
assert(max(values) <= edges(end) + sqrt(eps));
bin = discretize(values, edges);
assert(all(isfinite(bin)));
percentages = 100 * accumarray(bin, 1, ...
  [numel(edges) - 1 1], @sum, 0).' / numel(values);
end

function output = build_summary(state, catalog, parameters)
n = height(catalog);
sample = catalog.sample;
diameter_mm = catalog.diameter_mm;
cold_reduction_percent = catalog.cold_reduction_percent;
variant = catalog.variant;
input_path = catalog.input_path;
ebsd_step_um = reshape([state.step_um], n, 1);
kam_order = repmat(parameters.kam_order, n, 1);
kam_cutoff_deg = repmat(parameters.kam_cutoff_deg, n, 1);
burgers_vector_nm = repmat(parameters.burgers_vector_nm, n, 1);
valid_pixel_count = reshape([state.valid_pixel_count], n, 1);
kam_mean_deg = reshape([state.kam_mean_deg], n, 1);
kam_median_deg = reshape([state.kam_median_deg], n, 1);
kam_p90_deg = reshape([state.kam_p90_deg], n, 1);
gnd_mean_1e14_m2 = reshape([state.gnd_mean_1e14_m2], n, 1);
gnd_median_1e14_m2 = reshape([state.gnd_median_1e14_m2], n, 1);
gnd_p90_1e14_m2 = reshape([state.gnd_p90_1e14_m2], n, 1);
map_saturation_percent = reshape([state.map_saturation_percent], n, 1);
output = table(sample, diameter_mm, cold_reduction_percent, variant, ...
  input_path, ebsd_step_um, kam_order, kam_cutoff_deg, ...
  burgers_vector_nm, valid_pixel_count, kam_mean_deg, kam_median_deg, ...
  kam_p90_deg, gnd_mean_1e14_m2, gnd_median_1e14_m2, ...
  gnd_p90_1e14_m2, map_saturation_percent);
end

function output = build_distribution_table(state, catalog, edges)
nBins = numel(edges) - 1;
nRows = numel(state) * nBins;
sample = strings(nRows, 1);
diameter_mm = zeros(nRows, 1);
cold_reduction_percent = zeros(nRows, 1);
bin_lower_1e14_m2 = zeros(nRows, 1);
bin_upper_1e14_m2 = zeros(nRows, 1);
bin_center_1e14_m2 = zeros(nRows, 1);
pixel_frequency_percent = zeros(nRows, 1);
for stateIndex = 1:numel(state)
  rows = (stateIndex - 1) * nBins + (1:nBins);
  sample(rows) = repmat(catalog.sample(stateIndex), nBins, 1);
  diameter_mm(rows) = catalog.diameter_mm(stateIndex);
  cold_reduction_percent(rows) = ...
    catalog.cold_reduction_percent(stateIndex);
  bin_lower_1e14_m2(rows) = edges(1:end-1).';
  bin_upper_1e14_m2(rows) = edges(2:end).';
  bin_center_1e14_m2(rows) = ...
    (edges(1:end-1) + diff(edges) / 2).';
  pixel_frequency_percent(rows) = state(stateIndex).hist_percent(:);
end
output = table(sample, diameter_mm, cold_reduction_percent, ...
  bin_lower_1e14_m2, bin_upper_1e14_m2, bin_center_1e14_m2, ...
  pixel_frequency_percent);
end

function output = build_parameter_table(parameters, state)
parameter = [ ...
  "figure_basis"; "input_variant"; "grain_detection_deg"; ...
  "kam_order"; "kam_cutoff_deg"; ...
  "ebsd_step_um"; "burgers_vector_nm"; "gnd_equation"; ...
  "kam_angle_unit_in_equation"; "gnd_map_max_1e14_m2"; ...
  "gnd_hist_max_1e14_m2"; "gnd_hist_bin_width_1e14_m2"; ...
  "distribution_weighting"; "interpretive_limit"];
value = [ ...
  "layout and content adapted from Xia et al. MSEA 975 (2026) Fig. 14"; ...
  string(parameters.input_variant); ...
  string(parameters.grain_detection_deg); ...
  string(parameters.kam_order); ...
  string(parameters.kam_cutoff_deg); string(state(1).step_um); ...
  string(parameters.burgers_vector_nm); "rho_GND=2*KAM/(step*b)"; ...
  "radian"; string(parameters.gnd_map_max_1e14_m2); ...
  string(parameters.gnd_hist_max_1e14_m2); ...
  string(parameters.gnd_hist_bin_width_1e14_m2); ...
  "Ti-Hex pixel frequency"; ...
  "KAM-derived GND estimate; not total dislocation density"];
output = table(parameter, value);
end

function write_readme(outputRoot, parameters)
lines = [ ...
  "# Fig. 14参考形式：六直径KAM-GND图"; ""; ...
  "## 图件内容"; ""; ...
  "- 三行分别比较7-6.48 mm、6.02-5.6 mm和5.25-5 mm。"; ...
  "- 每行左侧为两张GND空间图，右侧为对应的上下两组像素频率分布。"; ...
  "- 六张空间图采用统一色标，不绘制等高线。"; ""; ...
  "## 计算方法"; ""; ...
  "- 输入为raw CTF；KAM采用一阶邻域和5 deg排除阈值。"; ...
  "- GND按rho_GND=2*KAM/(step*b)逐像素换算，其中KAM使用弧度。"; ...
  sprintf("- EBSD步长为0.5 um，Burgers矢量取%.3f nm。", ...
    parameters.burgers_vector_nm); ...
  "- 频率分布按有效Ti-Hex像素计数。"; ""; ...
  "## 解释边界"; ""; ...
  "该结果是基于KAM的几何必要位错密度估算，反映局部取向梯度，" + ...
  "不等同于包含统计存储位错在内的总位错密度。结果受步长、" + ...
  "邻域阶数、阈值、索引质量和数据清理方法影响。"; ""; ...
  "## 建议图注"; ""; ...
  "不同冷变形状态下Gr4B23271商业纯钛的KAM估算GND密度空间分布" + ...
  "及统计结果：（a-c）7和6.48 mm样品；（d-f）6.02和5.6 mm样品；" + ...
  "（g-i）5.25和5 mm样品。GND密度由一阶邻域KAM按" + ...
  "rho_GND=2*KAM/(step*b)估算，EBSD步长为0.5 um，" + ...
  "b=0.295 nm。空间图采用统一色标，统计分布按有效Ti-Hex像素计数。"; ...
  ""; "## 输出"; ""; ...
  "- `reference_style_fig14_gnd_six_diameters.*`：600 dpi PNG、TIFF和PDF。"; ...
  "- `reference_style_fig14_gnd_summary.csv`：六状态KAM和GND统计。"; ...
  "- `reference_style_fig14_gnd_distribution.csv`：完整分布数据。"; ...
  "- `reference_style_fig14_gnd_parameters.csv`：计算参数和解释边界。"];
readmePath = fullfile(outputRoot, "README.md");
fileId = fopen(readmePath, "w");
assert(fileId >= 0, "Could not open README for writing.");
cleanupFile = onCleanup(@() fclose(fileId));
fprintf(fileId, "%s\n", lines);
clear cleanupFile
end

function validate_outputs(state, summary, edges, paths)
assert(height(summary) == 6);
for stateIndex = 1:numel(state)
  assert(abs(sum(state(stateIndex).hist_percent) - 100) < 1e-8);
  assert(state(stateIndex).step_um > 0);
  assert(state(stateIndex).gnd_mean_1e14_m2 > 0);
end
assert(edges(1) == 0 && edges(end) == 12);
for pathIndex = 1:numel(paths)
  assert(isfile(paths(pathIndex)), "Missing output: %s", paths(pathIndex));
end
fprintf("REFERENCE_FIG14 validation passed: six states and distributions.\n");
end
