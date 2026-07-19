function [summary, sensitivity, distribution] = ...
  generate_grain_boundary_misorientation_distribution(scanRoot, outputDir)
%GENERATE_GRAIN_BOUNDARY_MISORIENTATION_DISTRIBUTION Analyze Ti-Hex boundaries.

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
detectionFloorsDeg = [0.5; 1; 2];
classificationAngleDeg = 15;
minPixel = 5;
minInnerComponentLengthUm = 1.0;
distributionEdgesDeg = 1:0.5:94;

sampleCount = numel(sampleNames);
floorCount = numel(detectionFloorsDeg);
summaryRows = cell(sampleCount, 1);
sensitivityRows = cell(sampleCount * floorCount, 1);
distributionBlocks = cell(sampleCount, 1);
sensitivityIndex = 0;

for sampleIndex = 1:sampleCount
  inputFile = fullfile(scanRoot, folderNames(sampleIndex), ...
    inputNames(sampleIndex));
  assert(isfile(inputFile), "CTF file not found: %s", inputFile);
  fprintf("IMPORT=%s\n", inputFile);

  ebsd = EBSD.load(inputFile, "convertEuler2SpatialReferenceFrame");
  indexedTiHexPixelCount = length(ebsd("Ti-Hex"));
  indexedTiHexAreaUm2 = indexedTiHexPixelCount * ...
    polyarea(ebsd.unitCell.x, ebsd.unitCell.y);
  assert(indexedTiHexPixelCount > 0 && ...
    isfinite(indexedTiHexAreaUm2) && indexedTiHexAreaUm2 > 0, ...
    "Ti-Hex indexed area is invalid for %s.", sampleNames(sampleIndex));

  for floorIndex = 1:floorCount
    floorDeg = detectionFloorsDeg(floorIndex);
    [grains, ebsd.grainId] = calcGrains(ebsd('indexed'), ...
      "threshold", [floorDeg classificationAngleDeg] * degree, ...
      "minPixel", minPixel);

    inner = grains.innerBoundary("Ti-Hex", "Ti-Hex");
    outer = grains.boundary("Ti-Hex", "Ti-Hex");
    innerTheta = double(angle(inner.misorientation) / degree);
    outerTheta = double(angle(outer.misorientation) / degree);
    innerTheta = innerTheta(:);
    outerTheta = outerTheta(:);
    innerLength = double(inner.segLength(:));
    outerLength = double(outer.segLength(:));

    if isempty(innerLength)
      keepInner = false(0,1);
    else
      [keepInner, ~] = component_length_mask(inner.componentId, ...
        innerLength, minInnerComponentLengthUm);
    end

    innerValid = isfinite(innerTheta) & isfinite(innerLength) & ...
      innerLength > 0 & innerTheta >= floorDeg - 1e-8 & innerTheta <= 94;
    outerValid = isfinite(outerTheta) & isfinite(outerLength) & ...
      outerLength > 0 & outerTheta >= floorDeg - 1e-8 & outerTheta <= 94;
    thetaDeg = [innerTheta(keepInner); outerTheta];
    segLengthUm = [innerLength(keepInner); outerLength];
    valid = isfinite(thetaDeg) & isfinite(segLengthUm) & ...
      segLengthUm > 0 & thetaDeg >= floorDeg - 1e-8 & thetaDeg <= 94;
    thetaDeg = thetaDeg(valid);
    segLengthUm = segLengthUm(valid);

    assert(~isempty(thetaDeg), ...
      "No retained Ti-Hex boundaries for %s at %.1f degrees.", ...
      sampleNames(sampleIndex), floorDeg);
    boundaryStats = calculate_boundary_metrics(thetaDeg, segLengthUm, ...
      classificationAngleDeg, indexedTiHexAreaUm2);

    rawInnerSegmentCount = nnz(innerValid);
    rawInnerBoundaryLengthUm = sum(innerLength(innerValid));
    retainedInnerMask = innerValid & keepInner;
    retainedInnerSegmentCount = nnz(retainedInnerMask);
    retainedInnerBoundaryLengthUm = sum(innerLength(retainedInnerMask));
    removedInnerSegmentCount = nnz(innerValid & ~keepInner);
    removedInnerBoundaryLengthUm = sum(innerLength(innerValid & ~keepInner));
    outerSegmentCount = nnz(outerValid);
    outerBoundaryLengthUm = sum(outerLength(outerValid));

    row = struct();
    row.sample = sampleNames(sampleIndex);
    row.folder = folderNames(sampleIndex);
    row.input_file = inputNames(sampleIndex);
    row.cold_reduction_percent = coldReductionPercent(sampleIndex);
    row.detection_floor_deg = floorDeg;
    row.classification_angle_deg = classificationAngleDeg;
    row.min_pixel = minPixel;
    row.min_inner_component_length_um = minInnerComponentLengthUm;
    row.grain_count = length(grains);
    row.indexed_ti_hex_pixel_count = indexedTiHexPixelCount;
    row.indexed_ti_hex_area_um2 = indexedTiHexAreaUm2;
    row.raw_inner_segment_count = rawInnerSegmentCount;
    row.raw_inner_boundary_length_um = rawInnerBoundaryLengthUm;
    row.retained_inner_segment_count = retainedInnerSegmentCount;
    row.retained_inner_boundary_length_um = retainedInnerBoundaryLengthUm;
    row.removed_inner_segment_count = removedInnerSegmentCount;
    row.removed_inner_boundary_length_um = removedInnerBoundaryLengthUm;
    row.outer_segment_count = outerSegmentCount;
    row.outer_boundary_length_um = outerBoundaryLengthUm;
    metricNames = fieldnames(boundaryStats);
    for metricIndex = 1:numel(metricNames)
      metricName = metricNames{metricIndex};
      row.(metricName) = boundaryStats.(metricName);
    end

    sensitivityIndex = sensitivityIndex + 1;
    sensitivityRows{sensitivityIndex} = struct2table(row);

    if floorDeg == 1
      intervalMasks = {
        thetaDeg >= 1 & thetaDeg < 2;
        thetaDeg >= 2 & thetaDeg < 5;
        thetaDeg >= 5 & thetaDeg < 10;
        thetaDeg >= 10 & thetaDeg < 15;
        thetaDeg >= 15 & thetaDeg <= 94
      };
      intervalLabels = ["1_2"; "2_5"; "5_10"; "10_15"; "15_94"];
      summaryRow = row;
      for intervalIndex = 1:numel(intervalMasks)
        intervalLength = sum(segLengthUm(intervalMasks{intervalIndex}));
        summaryRow.("length_" + intervalLabels(intervalIndex) + "_um") = ...
          intervalLength;
        summaryRow.("fraction_" + intervalLabels(intervalIndex)) = ...
          intervalLength / boundaryStats.total_retained_boundary_length_um;
      end
      summaryRows{sampleIndex} = struct2table(summaryRow);

      [weightedStats, sampleDistribution] = ...
        summarize_weighted_boundary_angles(thetaDeg, segLengthUm, ...
        distributionEdgesDeg);
      assert(abs(weightedStats.total_boundary_length_um - ...
        boundaryStats.total_retained_boundary_length_um) < 1e-8);
      distributionBlocks{sampleIndex} = addvars(sampleDistribution, ...
        repmat(sampleNames(sampleIndex), height(sampleDistribution), 1), ...
        repmat(folderNames(sampleIndex), height(sampleDistribution), 1), ...
        repmat(inputNames(sampleIndex), height(sampleDistribution), 1), ...
        repmat(coldReductionPercent(sampleIndex), ...
          height(sampleDistribution), 1), ...
        repmat(indexedTiHexAreaUm2, height(sampleDistribution), 1), ...
        'Before', 1, 'NewVariableNames', ...
        {'sample', 'folder', 'input_file', 'cold_reduction_percent', ...
        'indexed_ti_hex_area_um2'});
    end

    fprintf("SAMPLE=%s FLOOR_DEG=%.1f GRAINS=%d RETAINED_UM=%.6f " + ...
      "LAGB_FRACTION=%.6f REMOVED_INNER_UM=%.6f\n", ...
      sampleNames(sampleIndex), floorDeg, length(grains), ...
      boundaryStats.total_retained_boundary_length_um, ...
      boundaryStats.lagb_length_fraction, removedInnerBoundaryLengthUm);
    clear grains inner outer innerTheta outerTheta innerLength outerLength
  end
  clear ebsd
end

summary = vertcat(summaryRows{:});
sensitivity = vertcat(sensitivityRows{:});
distribution = vertcat(distributionBlocks{:});
validate_tables(summary, sensitivity, distribution, sampleNames, ...
  detectionFloorsDeg, distributionEdgesDeg);

summaryFile = fullfile(outputDir, ...
  "grain_boundary_misorientation_summary.csv");
sensitivityFile = fullfile(outputDir, ...
  "grain_boundary_detection_sensitivity.csv");
distributionFile = fullfile(outputDir, ...
  "grain_boundary_misorientation_distribution.csv");
distributionFigureFile = fullfile(outputDir, ...
  "grain_boundary_misorientation_distribution.png");
metricsFigureFile = fullfile(outputDir, ...
  "grain_boundary_misorientation_metrics.png");
writetable(summary, summaryFile);
writetable(sensitivity, sensitivityFile);
writetable(distribution, distributionFile);
plot_boundary_distributions(summary, distribution, distributionFigureFile);
plot_boundary_metrics(summary, metricsFigureFile);

fprintf("SUMMARY=%s\n", summaryFile);
fprintf("SENSITIVITY=%s\n", sensitivityFile);
fprintf("DISTRIBUTION=%s\n", distributionFile);
fprintf("DISTRIBUTION_FIGURE=%s\n", distributionFigureFile);
fprintf("METRICS_FIGURE=%s\n", metricsFigureFile);
end

function plot_boundary_distributions(summary, distribution, outputFile)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [100, 100, 1200, 900]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 2, 1, "Padding", "compact", ...
  "TileSpacing", "compact");
colors = lines(height(summary));
panelLimits = [1, 94; 1, 15];
panelTitles = ["Complete range (1-94 degree)"; ...
  "Low-angle detail (1-15 degree)"];

for panelIndex = 1:2
  axesHandle = nexttile(layout);
  hold(axesHandle, "on");
  for sampleIndex = 1:height(summary)
    rows = distribution.sample == summary.sample(sampleIndex);
    plot(axesHandle, distribution.bin_center_deg(rows), ...
      distribution.probability_density_per_degree(rows), ...
      "Color", colors(sampleIndex,:), "LineWidth", 1.7, ...
      "DisplayName", summary.sample(sampleIndex));
  end
  xlim(axesHandle, panelLimits(panelIndex,:));
  xlabel(axesHandle, ...
    "Ti-Hex grain-boundary misorientation (degree)");
  ylabel(axesHandle, ...
    "Length-weighted probability density (degree^{-1})");
  title(axesHandle, panelTitles(panelIndex));
  grid(axesHandle, "on");
  box(axesHandle, "on");
  text(axesHandle, 0.01, 0.92, sprintf("(%c)", 'a' + panelIndex - 1), ...
    "Units", "normalized", "FontWeight", "bold");
  if panelIndex == 1
    legend(axesHandle, "Location", "eastoutside");
  end
end

export_clean_png(figureHandle, outputFile);
clear cleanupFigure
end

function plot_boundary_metrics(summary, outputFile)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [100, 100, 1200, 800]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 2, 1, "Padding", "compact", ...
  "TileSpacing", "compact");
colors = lines(height(summary));
xValues = summary.cold_reduction_percent;
yValues = [100 * summary.lagb_length_fraction, ...
  summary.lagb_length_density_um_per_um2];
yLabels = ["LAGB length fraction (%)"; ...
  "LAGB length density (\mum \mum^{-2})"];
panelTitles = ["Low-angle grain-boundary length fraction"; ...
  "Low-angle grain-boundary length density"];

for panelIndex = 1:2
  axesHandle = nexttile(layout);
  hold(axesHandle, "on");
  plot(axesHandle, xValues, yValues(:,panelIndex), "-", ...
    "Color", [0.35, 0.35, 0.35], "LineWidth", 1.2, ...
    "HandleVisibility", "off");
  for sampleIndex = 1:height(summary)
    plot(axesHandle, xValues(sampleIndex), ...
      yValues(sampleIndex,panelIndex), "o", ...
      "Color", colors(sampleIndex,:), ...
      "MarkerFaceColor", colors(sampleIndex,:), "MarkerSize", 7, ...
      "HandleVisibility", "off");
    text(axesHandle, xValues(sampleIndex), ...
      yValues(sampleIndex,panelIndex), "  " + summary.sample(sampleIndex), ...
      "HorizontalAlignment", "left", "VerticalAlignment", "bottom");
  end
  xlim(axesHandle, [-2, 55]);
  xlabel(axesHandle, "Cold reduction (%)");
  ylabel(axesHandle, yLabels(panelIndex));
  title(axesHandle, panelTitles(panelIndex));
  grid(axesHandle, "on");
  box(axesHandle, "on");
  text(axesHandle, 0.01, 0.92, sprintf("(%c)", 'a' + panelIndex - 1), ...
    "Units", "normalized", "FontWeight", "bold");
end

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

function stats = calculate_boundary_metrics(thetaDeg, segLengthUm, ...
  classificationAngleDeg, indexedAreaUm2)
totalLengthUm = sum(segLengthUm);
lagbMask = thetaDeg < classificationAngleDeg;
hagbMask = thetaDeg >= classificationAngleDeg;
lagbLengthUm = sum(segLengthUm(lagbMask));
hagbLengthUm = sum(segLengthUm(hagbMask));
[sortedTheta, order] = sort(thetaDeg);
sortedLength = segLengthUm(order);
weightedMedianIndex = find(cumsum(sortedLength) >= 0.5 * totalLengthUm, ...
  1);

stats.retained_segment_count = numel(thetaDeg);
stats.total_retained_boundary_length_um = totalLengthUm;
stats.lagb_segment_count = nnz(lagbMask);
stats.hagb_segment_count = nnz(hagbMask);
stats.lagb_length_um = lagbLengthUm;
stats.hagb_length_um = hagbLengthUm;
stats.lagb_length_fraction = lagbLengthUm / totalLengthUm;
stats.hagb_length_fraction = hagbLengthUm / totalLengthUm;
stats.lagb_length_density_um_per_um2 = lagbLengthUm / indexedAreaUm2;
stats.hagb_length_density_um_per_um2 = hagbLengthUm / indexedAreaUm2;
stats.weighted_mean_angle_deg = sum(thetaDeg .* segLengthUm) / totalLengthUm;
stats.weighted_median_angle_deg = sortedTheta(weightedMedianIndex);
stats.min_angle_deg = min(thetaDeg);
stats.max_angle_deg = max(thetaDeg);
end

function validate_tables(summary, sensitivity, distribution, sampleNames, ...
  detectionFloorsDeg, distributionEdgesDeg)
expectedDistributionRows = numel(sampleNames) * ...
  (numel(distributionEdgesDeg) - 1);
assert(height(summary) == 6, "Expected six primary summary rows.");
assert(height(sensitivity) == 18, "Expected eighteen sensitivity rows.");
assert(height(distribution) == expectedDistributionRows, ...
  "Expected 1116 distribution rows.");
assert(isequal(summary.sample, sampleNames), ...
  "Primary samples are missing or out of order.");
assert(isequal(unique(sensitivity.detection_floor_deg), ...
  detectionFloorsDeg), "Detection floors are missing.");

assert(all(isfinite(sensitivity.indexed_ti_hex_area_um2) & ...
  sensitivity.indexed_ti_hex_area_um2 > 0), ...
  "Indexed Ti-Hex areas must be positive finite values.");
assert(all(isfinite(sensitivity.total_retained_boundary_length_um) & ...
  sensitivity.total_retained_boundary_length_um > 0), ...
  "Retained boundary lengths must be positive finite values.");
assert(all(abs(sensitivity.lagb_length_um + sensitivity.hagb_length_um - ...
  sensitivity.total_retained_boundary_length_um) < 1e-8), ...
  "LAGB and HAGB lengths do not conserve retained length.");
assert(all(abs(sensitivity.lagb_length_fraction + ...
  sensitivity.hagb_length_fraction - 1) < 1e-10), ...
  "LAGB and HAGB fractions do not sum to one.");
assert(all(sensitivity.lagb_length_fraction >= 0 & ...
  sensitivity.lagb_length_fraction <= 1 & ...
  sensitivity.hagb_length_fraction >= 0 & ...
  sensitivity.hagb_length_fraction <= 1), ...
  "Boundary fractions fall outside [0,1].");
assert(all(isfinite(sensitivity.weighted_mean_angle_deg) & ...
  isfinite(sensitivity.weighted_median_angle_deg) & ...
  isfinite(sensitivity.min_angle_deg) & ...
  isfinite(sensitivity.max_angle_deg)), ...
  "Boundary angle summaries must be finite.");
assert(all(abs(sensitivity.raw_inner_boundary_length_um - ...
  sensitivity.retained_inner_boundary_length_um - ...
  sensitivity.removed_inner_boundary_length_um) < 1e-8), ...
  "Inner-boundary continuity-filter audit does not conserve length.");
assert(all(abs(sensitivity.retained_inner_boundary_length_um + ...
  sensitivity.outer_boundary_length_um - ...
  sensitivity.total_retained_boundary_length_um) < 1e-8), ...
  "Retained inner and outer lengths do not conserve total length.");

intervalLength = summary.length_1_2_um + summary.length_2_5_um + ...
  summary.length_5_10_um + summary.length_10_15_um + ...
  summary.length_15_94_um;
assert(all(abs(intervalLength - ...
  summary.total_retained_boundary_length_um) < 1e-8), ...
  "Primary interval lengths do not conserve retained length.");
intervalFraction = summary.fraction_1_2 + summary.fraction_2_5 + ...
  summary.fraction_5_10 + summary.fraction_10_15 + ...
  summary.fraction_15_94;
assert(all(abs(intervalFraction - 1) < 1e-10), ...
  "Primary interval fractions do not sum to one.");

for sampleIndex = 1:numel(sampleNames)
  sensitivityMask = sensitivity.sample == sampleNames(sampleIndex);
  assert(nnz(sensitivityMask) == numel(detectionFloorsDeg), ...
    "Expected three sensitivity rows for %s.", sampleNames(sampleIndex));
  assert(isequal(sort(sensitivity.detection_floor_deg(sensitivityMask)), ...
    detectionFloorsDeg), "Detection floors are incomplete for %s.", ...
    sampleNames(sampleIndex));
  distributionMask = distribution.sample == sampleNames(sampleIndex);
  assert(nnz(distributionMask) == numel(distributionEdgesDeg) - 1, ...
    "Unexpected distribution bin count for %s.", sampleNames(sampleIndex));
  assert(abs(sum(distribution.length_fraction(distributionMask)) - 1) ...
    < 1e-10, "Length fractions do not sum to one for %s.", ...
    sampleNames(sampleIndex));
  assert(abs(sum(distribution.probability_density_per_degree( ...
    distributionMask) .* distribution.bin_width_deg(distributionMask)) ...
    - 1) < 1e-10, "Probability density does not integrate to one for %s.", ...
    sampleNames(sampleIndex));
end
end
