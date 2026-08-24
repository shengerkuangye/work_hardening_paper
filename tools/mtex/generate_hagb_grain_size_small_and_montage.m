function generate_hagb_grain_size_small_and_montage(sourceDir, outputRoot)
%GENERATE_HAGB_GRAIN_SIZE_SMALL_AND_MONTAGE Plot conventional grain ECD.
% Uses the previously calculated raw-EBSD, 15 degree HAGB grain
% reconstruction (minPixel = 5). Exports six standalone 600-dpi PNG files
% and one common-axis 2-by-3 montage in PNG, TIFF, and PDF formats.

arguments
  sourceDir (1,1) string
  outputRoot (1,1) string
end

summaryPath = fullfile(sourceDir, "hagb_grain_size_summary.csv");
histogramPath = fullfile(sourceDir, "hagb_grain_size_histograms.csv");
assert(isfile(summaryPath), "Missing HAGB grain-size summary: %s", summaryPath);
assert(isfile(histogramPath), "Missing HAGB histogram data: %s", histogramPath);

summary = readtable(summaryPath, "TextType", "string");
histograms = readtable(histogramPath, "TextType", "string");
[~, order] = sort(summary.cold_reduction_percent);
summary = summary(order, :);
assert(height(summary) == 6, "Expected six deformation states.");
assert(all(summary.hagb_threshold_deg == 15));
assert(all(summary.min_grain_pixels == 5));

smallDir = fullfile(outputRoot, "07_hagb_grain_size");
montageDir = fullfile(outputRoot, "montages");
if ~isfolder(smallDir), mkdir(smallDir); end
if ~isfolder(montageDir), mkdir(montageDir); end

colors = [ ...
  0.40 0.40 0.40
  0.00 0.45 0.70
  0.34 0.71 0.91
  0.00 0.62 0.45
  0.90 0.62 0.00
  0.84 0.37 0.00];

state = repmat(empty_state(), 6, 1);
for stateIndex = 1:6
  rows = histograms(histograms.sample == summary.sample(stateIndex), :);
  [~, rowOrder] = sort(rows.bin_center_um);
  rows = rows(rowOrder, :);
  state(stateIndex).edges = [rows.bin_lower_um(1); rows.bin_upper_um];
  state(stateIndex).centers = rows.bin_center_um(:).';
  state(stateIndex).frequency = rows.number_frequency_percent(:).';
  state(stateIndex).mean_um = summary.ecd_number_mean_um(stateIndex);
  state(stateIndex).median_um = summary.ecd_number_median_um(stateIndex);
  state(stateIndex).mu = summary.number_lognormal_mu(stateIndex);
  state(stateIndex).sigma = summary.number_lognormal_sigma(stateIndex);
  state(stateIndex).grain_count = summary.grain_count(stateIndex);
  assert(abs(sum(state(stateIndex).frequency) - 100) < 1e-8);
end

assert(all(arrayfun(@(s) isequal(s.edges, state(1).edges), state)), ...
  "All states must use identical HAGB ECD bins.");
commonEdges = state(1).edges;
commonYMax = calculate_ymax(state, commonEdges);

for stateIndex = 1:6
  figureHandle = figure("Visible", "off", "Color", "white", ...
    "Units", "inches", "Position", [0.25 0.25 4.2 3.2]);
  cleanupFigure = onCleanup(@() close(figureHandle));
  axesHandle = axes(figureHandle, "Units", "normalized", ...
    "Position", [0.17 0.18 0.78 0.69]);
  localYMax = calculate_ymax(state(stateIndex), commonEdges);
  draw_panel(axesHandle, state(stateIndex), colors(stateIndex, :), ...
    commonEdges, localYMax, true);
  title(axesHandle, sprintf("%.2f%% cold reduction | %.2f mm", ...
    summary.cold_reduction_percent(stateIndex), ...
    summary.diameter_mm(stateIndex)), "FontName", "Arial", ...
    "FontSize", 9, "FontWeight", "bold");
  tag = state_tag(summary(stateIndex, :));
  exportgraphics(figureHandle, fullfile(smallDir, ...
    tag + "_hagb_grain_ecd.png"), "Resolution", 600, ...
    "BackgroundColor", "white");
  clear cleanupFigure
end

figureHandle = figure("Visible", "off", "Color", "white", ...
  "Units", "inches", "Position", [0.25 0.25 11.2 7.2]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 2, 3, "Padding", "compact", ...
  "TileSpacing", "compact");
for stateIndex = 1:6
  axesHandle = nexttile(layout);
  draw_panel(axesHandle, state(stateIndex), colors(stateIndex, :), ...
    commonEdges, commonYMax, stateIndex == 1);
  title(axesHandle, sprintf("(%c) %.2f%% | %.2f mm", ...
    char('a' + stateIndex - 1), ...
    summary.cold_reduction_percent(stateIndex), ...
    summary.diameter_mm(stateIndex)), "FontName", "Arial", ...
    "FontSize", 8, "FontWeight", "bold");
end
title(layout, "Conventional grain ECD | 15 deg HAGB; number frequency; minPixel = 5", ...
  "FontName", "Arial", "FontSize", 10, "FontWeight", "bold");
basePath = fullfile(montageDir, "hagb_grain_size_six_state_montage");
exportgraphics(figureHandle, basePath + ".png", "Resolution", 600, ...
  "BackgroundColor", "white");
exportgraphics(figureHandle, basePath + ".tif", "Resolution", 600, ...
  "BackgroundColor", "white");
exportgraphics(figureHandle, basePath + ".pdf", ...
  "ContentType", "vector", "BackgroundColor", "white");
clear cleanupFigure

writetable(summary, fullfile(smallDir, "hagb_grain_size_summary_used.csv"));
assert(numel(dir(fullfile(smallDir, "*_hagb_grain_ecd.png"))) == 6);
fprintf("HAGB_GRAIN_SIZE_GALLERY validation passed: 6 panels + 1 montage.\n");
end

function result = empty_state()
result = struct("edges", [], "centers", [], "frequency", [], ...
  "mean_um", NaN, "median_um", NaN, "mu", NaN, "sigma", NaN, ...
  "grain_count", 0);
end

function draw_panel(axesHandle, state, color, edges, yMaximum, showLegend)
barHandle = bar(axesHandle, state.centers, state.frequency, 0.84, ...
  "FaceColor", color, "FaceAlpha", 0.74, ...
  "EdgeColor", [0.16 0.16 0.16], "LineWidth", 0.35);
hold(axesHandle, "on");
fineX = linspace(max(0.05, edges(1) + 0.05), edges(end), 600);
fineY = 100 * median(diff(edges)) * lognormal_pdf(fineX, ...
  state.mu, state.sigma);
fitHandle = plot(axesHandle, fineX, fineY, ...
  "Color", [0.08 0.08 0.08], "LineWidth", 1.35);
if showLegend
  legend(axesHandle, [barHandle fitHandle], ...
    ["Histogram", "Descriptive lognormal fit"], ...
    "Location", "northeast", "FontName", "Arial", ...
    "FontSize", 6.5, "Box", "off");
  meanTextY = 0.70;
else
  meanTextY = 0.92;
end
text(axesHandle, 0.97, meanTextY, sprintf("Mean ECD = %.2f um", ...
  state.mean_um), "Units", "normalized", ...
  "HorizontalAlignment", "right", "VerticalAlignment", "top", ...
  "FontName", "Arial", "FontSize", 7.5, "FontWeight", "bold", ...
  "Color", color);
xlim(axesHandle, [edges(1) edges(end)]);
ylim(axesHandle, [0 yMaximum]);
xlabel(axesHandle, "15 deg HAGB grain ECD (um)", ...
  "FontName", "Arial", "FontSize", 8);
ylabel(axesHandle, "Number frequency (%)", ...
  "FontName", "Arial", "FontSize", 8);
set(axesHandle, "FontName", "Arial", "FontSize", 7, ...
  "TickDir", "out", "Box", "on", "LineWidth", 0.6, "Layer", "top");
end

function value = calculate_ymax(state, edges)
fineX = linspace(max(0.05, edges(1) + 0.05), edges(end), 600);
peak = 0;
for stateIndex = 1:numel(state)
  fitted = 100 * median(diff(edges)) * lognormal_pdf(fineX, ...
    state(stateIndex).mu, state(stateIndex).sigma);
  peak = max([peak, state(stateIndex).frequency, fitted]);
end
value = max(5, 5 * ceil(1.12 * peak / 5));
end

function density = lognormal_pdf(values, mu, sigma)
density = exp(-0.5 * ((log(values) - mu) / sigma).^2) ./ ...
  (values * sigma * sqrt(2 * pi));
end

function tag = state_tag(summaryRow)
diameter = replace(sprintf("%.2f", summaryRow.diameter_mm), ".", "p");
reduction = replace(sprintf("%.2f", ...
  summaryRow.cold_reduction_percent), ".", "p");
tag = diameter + "mm_" + reduction + "pct";
end
