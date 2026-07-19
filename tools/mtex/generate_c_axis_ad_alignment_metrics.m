function [summary, distribution] = generate_c_axis_ad_alignment_metrics( ...
  scanRoot, outputDir)
%GENERATE_C_AXIS_AD_ALIGNMENT_METRICS Analyze alpha-Ti c-axis angles to AD.

arguments
  scanRoot (1,1) string
  outputDir (1,1) string
end

assert(isfolder(scanRoot), "EBSD scan folder not found: %s", scanRoot);
assert(~isempty(which("EBSD")), "MTEX is not loaded in the current MATLAB session.");

if ~isfolder(outputDir)
  mkdir(outputDir);
end

folderNames = ["d7"; "d6_48"; "d6_02"; "d5_6"; "d5_25"; "d5"];
inputNames = [
  "ebsd_sample_7_map_15.ctf";
  "ebsd_sample_648_map_13.ctf";
  "ebsd_sample_602_map_11.ctf";
  "ebsd_sample_56_map_9.ctf";
  "ebsd_sample_525_map_7.ctf";
  "ebsd_sample_5_map_3.ctf"
];
sampleNames = ["7d"; "6.48d"; "6.02d"; "5.6d"; "5.25d"; "5d"];
coldReductionPercent = [0; 14.31; 26.04; 36.00; 43.75; 48.98];
binEdgesDeg = 0:2:90;

randomMeanAngleDeg = 180 / pi;
randomMedianAngleDeg = 60;
randomAlignmentFactor = 0;
randomFraction15 = 1 - cosd(15);
randomFraction30 = 1 - cosd(30);
randomFraction45 = 1 - cosd(45);

n = numel(sampleNames);
summary = table();
distribution = table();

for i = 1:n
  inputFile = fullfile(scanRoot, folderNames(i), inputNames(i));
  assert(isfile(inputFile), "CTF file not found: %s", inputFile);

  fprintf("IMPORT=%s\n", inputFile);
  ebsd = EBSD.load(inputFile, "convertEuler2SpatialReferenceFrame");
  alpha = ebsd("Ti-Hex");
  assert(~isempty(alpha), "Ti-Hex phase missing in %s", inputFile);

  cAxis = Miller(0, 0, 0, 1, alpha.CS);
  cDirections = alpha.orientations * cAxis;
  thetaDeg = c_axis_angles_to_ad(cDirections);
  [sampleStats, sampleDistribution] = compute_c_axis_ad_statistics( ...
    thetaDeg, length(ebsd), binEdgesDeg);

  summaryRow = struct2table(sampleStats);
  summaryRow = addvars(summaryRow, sampleNames(i), folderNames(i), ...
    inputNames(i), coldReductionPercent(i), 'Before', 1, ...
    'NewVariableNames', {'sample', 'folder', 'input_file', ...
    'cold_reduction_percent'});
  summaryRow = addvars(summaryRow, randomMeanAngleDeg, ...
    randomMedianAngleDeg, randomAlignmentFactor, randomFraction15, ...
    randomFraction30, randomFraction45, ...
    'NewVariableNames', {'random_mean_angle_deg', ...
    'random_median_angle_deg', 'random_alignment_factor', ...
    'random_fraction_within_15deg', 'random_fraction_within_30deg', ...
    'random_fraction_within_45deg'});

  sampleDistribution = addvars(sampleDistribution, ...
    repmat(sampleNames(i), height(sampleDistribution), 1), ...
    repmat(coldReductionPercent(i), height(sampleDistribution), 1), ...
    'Before', 1, 'NewVariableNames', ...
    {'sample', 'cold_reduction_percent'});

  if i == 1
    summary = summaryRow;
    distribution = sampleDistribution;
  else
    summary = [summary; summaryRow]; %#ok<AGROW>
    distribution = [distribution; sampleDistribution]; %#ok<AGROW>
  end

  fprintf("SAMPLE=%s ALPHA_PIXELS=%d MEAN_ANGLE_DEG=%.4f F_AD=%.6f\n", ...
    sampleNames(i), sampleStats.alpha_ti_pixel_count, ...
    sampleStats.mean_angle_deg, sampleStats.alignment_factor);

  clear ebsd alpha cAxis cDirections thetaDeg
end

assert(height(summary) == n, "Expected six summary rows.");
assert(height(distribution) == n * (numel(binEdgesDeg) - 1), ...
  "Unexpected angle-distribution row count.");
assert(all(isfinite(summary{:, 5:end}), "all"), ...
  "Summary contains nonfinite numeric values.");
assert(all(summary.alignment_factor >= -0.5 & ...
  summary.alignment_factor <= 1), "Alignment factor is outside its valid range.");

for i = 1:n
  rows = distribution.sample == sampleNames(i);
  assert(sum(distribution.pixel_count(rows)) == ...
    summary.alpha_ti_pixel_count(i), "Binned pixel count mismatch for %s.", ...
    sampleNames(i));
  assert(abs(sum(distribution.probability(rows)) - 1) < 1e-10, ...
    "Probability sum mismatch for %s.", sampleNames(i));
  assert(abs(sum(distribution.probability_density_per_degree(rows) .* ...
    distribution.bin_width_deg(rows)) - 1) < 1e-10, ...
    "Probability-density integral mismatch for %s.", sampleNames(i));
end

summaryFile = fullfile(outputDir, "c_axis_ad_alignment_summary.csv");
distributionFile = fullfile(outputDir, "c_axis_ad_angle_distribution.csv");
writetable(summary, summaryFile);
writetable(distribution, distributionFile);

metricsFigureFile = fullfile(outputDir, "c_axis_ad_alignment_metrics.png");
distributionFigureFile = fullfile(outputDir, ...
  "c_axis_ad_angle_distribution.png");
plot_alignment_metrics(summary, metricsFigureFile);
plot_angle_distribution(summary, distribution, distributionFigureFile);

fprintf("SUMMARY=%s\n", summaryFile);
fprintf("DISTRIBUTION=%s\n", distributionFile);
fprintf("METRICS_FIGURE=%s\n", metricsFigureFile);
fprintf("DISTRIBUTION_FIGURE=%s\n", distributionFigureFile);
end

function plot_alignment_metrics(summary, outputFile)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [100, 100, 1200, 1000]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 3, 1, "Padding", "compact", ...
  "TileSpacing", "compact");
x = summary.cold_reduction_percent;

axesHandle = nexttile(layout);
plot(axesHandle, x, summary.alignment_factor, "-o", ...
  "LineWidth", 1.8, "MarkerFaceColor", "auto");
hold(axesHandle, "on");
yline(axesHandle, 0, "--", "Random", "LineWidth", 1.2);
ylim(axesHandle, [-0.5, 0.05]);
ylabel(axesHandle, "F_{AD}");
grid(axesHandle, "on");
text(axesHandle, 0.01, 0.90, "(a)", "Units", "normalized", ...
  "FontWeight", "bold");

axesHandle = nexttile(layout);
plot(axesHandle, x, summary.mean_angle_deg, "-o", ...
  "LineWidth", 1.8, "MarkerFaceColor", "auto", ...
  "DisplayName", "Mean");
hold(axesHandle, "on");
plot(axesHandle, x, summary.median_angle_deg, "-s", ...
  "LineWidth", 1.8, "MarkerFaceColor", "auto", ...
  "DisplayName", "Median");
yline(axesHandle, 180/pi, "--", "Random mean", "LineWidth", 1.2, ...
  "HandleVisibility", "off");
yline(axesHandle, 60, ":", "Random median", "LineWidth", 1.2, ...
  "HandleVisibility", "off");
ylabel(axesHandle, "Angle to AD (degree)");
legend(axesHandle, "Location", "best");
grid(axesHandle, "on");
text(axesHandle, 0.01, 0.90, "(b)", "Units", "normalized", ...
  "FontWeight", "bold");

axesHandle = nexttile(layout);
plot(axesHandle, x, 100 * summary.fraction_within_15deg, "-o", ...
  x, 100 * summary.fraction_within_30deg, "-s", ...
  x, 100 * summary.fraction_within_45deg, "-^", ...
  "LineWidth", 1.8, "MarkerFaceColor", "auto");
ylabel(axesHandle, "Fraction (%)");
xticks(axesHandle, x);
xticklabels(axesHandle, compose("%g%% (%s)", x, summary.sample));
xtickangle(axesHandle, 20);
xlabel(axesHandle, "Cold reduction / sample");
legend(axesHandle, {"<=15 degree", "<=30 degree", "<=45 degree"}, ...
  "Location", "best");
grid(axesHandle, "on");
text(axesHandle, 0.01, 0.90, "(c)", "Units", "normalized", ...
  "FontWeight", "bold");

title(layout, "Alpha-Ti c-axis alignment relative to bar AD");
export_clean_png(figureHandle, outputFile);
clear cleanupFigure
end

function plot_angle_distribution(summary, distribution, outputFile)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [100, 100, 1200, 800]);
cleanupFigure = onCleanup(@() close(figureHandle));
axesHandle = axes(figureHandle);
hold(axesHandle, "on");
colors = lines(height(summary));

for i = 1:height(summary)
  rows = distribution.sample == summary.sample(i);
  plot(axesHandle, distribution.bin_center_deg(rows), ...
    distribution.probability_density_per_degree(rows), ...
    "LineWidth", 1.8, "Color", colors(i,:), ...
    "DisplayName", summary.sample(i));
end

thetaReference = linspace(0, 90, 361)';
randomDensityPerDegree = (pi/180) * sind(thetaReference);
plot(axesHandle, thetaReference, randomDensityPerDegree, "k--", ...
  "LineWidth", 1.8, "DisplayName", "Random");
xlim(axesHandle, [0, 90]);
xlabel(axesHandle, "c-axis angle to AD (degree)");
ylabel(axesHandle, "Probability density (degree^{-1})");
title(axesHandle, "Alpha-Ti c-axis angle distribution relative to bar AD");
legend(axesHandle, "Location", "eastoutside");
grid(axesHandle, "on");
export_clean_png(figureHandle, outputFile);
clear cleanupFigure
end

function export_clean_png(figureHandle, outputFile)
axesHandles = findall(figureHandle, "Type", "axes");
for axesHandle = reshape(axesHandles, 1, [])
  if isprop(axesHandle, "Toolbar")
    axesHandle.Toolbar = [];
  end
end
drawnow;
exportgraphics(figureHandle, char(outputFile), "Resolution", 300, ...
  "BackgroundColor", "white");
renderedImage = imread(outputFile);
imwrite(renderedImage, outputFile, "png");
end
