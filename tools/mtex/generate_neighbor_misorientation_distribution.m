function [summary, distribution] = ...
  generate_neighbor_misorientation_distribution(scanRoot, outputDir)
%GENERATE_NEIGHBOR_MISORIENTATION_DISTRIBUTION Compare Ti-Hex neighbors.

arguments
  scanRoot (1,1) string
  outputDir (1,1) string
end

assert(isfolder(scanRoot), "EBSD scan folder not found: %s", scanRoot);
assert(~isempty(which("EBSD")), ...
  "MTEX is not loaded in the current MATLAB session.");
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
binEdgesDeg = 0:0.1:94;

n = numel(sampleNames);
summaryRows = cell(n + 1, 1);
distributionBlocks = cell(n + 1, 1);
angleVectors = cell(n, 1);

for i = 1:n
  inputFile = fullfile(scanRoot, folderNames(i), inputNames(i));
  assert(isfile(inputFile), "CTF file not found: %s", inputFile);
  fprintf("IMPORT=%s\n", inputFile);

  ebsd = EBSD.load(inputFile, "convertEuler2SpatialReferenceFrame");
  ebsdGrid = ebsd.gridify;
  alpha = ebsdGrid("Ti-Hex");
  assert(~isempty(alpha), "Ti-Hex phase missing in %s", inputFile);
  alphaPhaseId = unique(alpha.phaseId);
  assert(isscalar(alphaPhaseId), ...
    "Expected one Ti-Hex phase identifier in %s", inputFile);

  phaseMask = reshape(ebsdGrid.phaseId, size(ebsdGrid)) == alphaPhaseId;
  [index1, index2] = four_neighbor_pairs(phaseMask);
  orientation1 = ebsdGrid(index1).orientations;
  orientation2 = ebsdGrid(index2).orientations;
  thetaDeg = double(angle(orientation1, orientation2) / degree);
  thetaDeg = thetaDeg(:);
  [sampleStats, sampleDistribution] = ...
    summarize_misorientation_angles(thetaDeg, binEdgesDeg);

  summaryRow = struct2table(sampleStats);
  summaryRows{i} = addvars(summaryRow, sampleNames(i), folderNames(i), ...
    inputNames(i), coldReductionPercent(i), 'Before', 1, ...
    'NewVariableNames', {'sample', 'folder', 'input_file', ...
    'cold_reduction_percent'});
  distributionBlocks{i} = addvars(sampleDistribution, ...
    repmat(sampleNames(i), height(sampleDistribution), 1), ...
    repmat(coldReductionPercent(i), height(sampleDistribution), 1), ...
    'Before', 1, 'NewVariableNames', ...
    {'sample', 'cold_reduction_percent'});
  angleVectors{i} = thetaDeg;

  fprintf("SAMPLE=%s PAIRS=%d MEAN_DEG=%.4f MEDIAN_DEG=%.4f MODE_BIN_DEG=%.2f\n", ...
    sampleNames(i), sampleStats.pair_count, ...
    sampleStats.mean_angle_deg, sampleStats.median_angle_deg, ...
    sampleStats.mode_bin_center_deg);

  clear ebsd ebsdGrid alpha phaseMask index1 index2
  clear orientation1 orientation2 thetaDeg sampleDistribution
end

pooledAngles = vertcat(angleVectors{:});
[pooledStats, pooledDistribution] = ...
  summarize_misorientation_angles(pooledAngles, binEdgesDeg);
pooledSummaryRow = struct2table(pooledStats);
summaryRows{n + 1} = addvars(pooledSummaryRow, "pooled", "", "", NaN, ...
  'Before', 1, 'NewVariableNames', {'sample', 'folder', 'input_file', ...
  'cold_reduction_percent'});
distributionBlocks{n + 1} = addvars(pooledDistribution, ...
  repmat("pooled", height(pooledDistribution), 1), ...
  NaN(height(pooledDistribution), 1), 'Before', 1, ...
  'NewVariableNames', {'sample', 'cold_reduction_percent'});

summary = vertcat(summaryRows{:});
distribution = vertcat(distributionBlocks{:});
validate_results(summary, distribution, sampleNames, binEdgesDeg);

summaryFile = fullfile(outputDir, "neighbor_misorientation_summary.csv");
distributionFile = fullfile(outputDir, ...
  "neighbor_misorientation_distribution.csv");
figureFile = fullfile(outputDir, ...
  "neighbor_misorientation_distribution.png");
writetable(summary, summaryFile);
writetable(distribution, distributionFile);
plot_distributions(summary, distribution, figureFile);

fprintf("SUMMARY=%s\n", summaryFile);
fprintf("DISTRIBUTION=%s\n", distributionFile);
fprintf("FIGURE=%s\n", figureFile);
end

function validate_results(summary, distribution, sampleNames, binEdgesDeg)
binCount = numel(binEdgesDeg) - 1;
assert(height(summary) == 7, "Expected six samples and one pooled row.");
assert(height(distribution) == 7 * binCount, ...
  "Unexpected distribution row count.");
assert(all(summary.pair_count(1:6) > 0 & ...
  summary.pair_count(1:6) <= 718800), "Invalid sample pair count.");
assert(summary.pair_count(7) == sum(summary.pair_count(1:6)), ...
  "Pooled pair count does not match the sample sum.");

for i = 1:numel(sampleNames)
  rows = distribution.sample == sampleNames(i);
  assert(nnz(rows) == binCount, "Unexpected bin count for %s.", ...
    sampleNames(i));
  assert(sum(distribution.pair_count(rows)) == summary.pair_count(i), ...
    "Binned pair count mismatch for %s.", sampleNames(i));
  assert(abs(sum(distribution.probability(rows)) - 1) < 1e-10, ...
    "Probability sum mismatch for %s.", sampleNames(i));
end

sampleCounts = reshape(distribution.pair_count(1:6*binCount), ...
  binCount, 6);
pooledRows = distribution.sample == "pooled";
assert(isequal(distribution.pair_count(pooledRows), ...
  sum(sampleCounts, 2)), "Pooled bin counts do not match sample sums.");
assert(abs(sum(distribution.probability(pooledRows)) - 1) < 1e-10, ...
  "Pooled probability sum mismatch.");
end

function plot_distributions(summary, distribution, outputFile)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [100, 100, 1200, 900]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 2, 1, "Padding", "compact", ...
  "TileSpacing", "compact");
colors = lines(6);

for panel = 1:2
  axesHandle = nexttile(layout);
  hold(axesHandle, "on");
  for i = 1:6
    rows = distribution.sample == summary.sample(i);
    plot(axesHandle, distribution.bin_center_deg(rows), ...
      distribution.probability_density_per_degree(rows), ...
      "Color", colors(i,:), "LineWidth", 1.5, ...
      "DisplayName", summary.sample(i));
  end
  pooledRows = distribution.sample == "pooled";
  plot(axesHandle, distribution.bin_center_deg(pooledRows), ...
    distribution.probability_density_per_degree(pooledRows), ...
    "k--", "LineWidth", 2.2, "DisplayName", "Pooled");
  ylabel(axesHandle, "Probability density (degree^{-1})");
  grid(axesHandle, "on");
  if panel == 1
    xlim(axesHandle, [0, 94]);
    title(axesHandle, "Full range (0-94 degree)");
    legend(axesHandle, "Location", "eastoutside");
    text(axesHandle, 0.01, 0.91, "(a)", "Units", "normalized", ...
      "FontWeight", "bold");
  else
    xlim(axesHandle, [0, 10]);
    title(axesHandle, "Low-angle detail (0-10 degree)");
    xlabel(axesHandle, "Ti-Hex four-neighbor misorientation (degree)");
    text(axesHandle, 0.01, 0.91, "(b)", "Units", "normalized", ...
      "FontWeight", "bold");
  end
end
title(layout, ["Symmetry-reduced Ti-Hex nearest-neighbor " ...
  "misorientation distribution"]);
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
