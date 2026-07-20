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

artifactNames = [
  "grain_boundary_detection_sensitivity.csv";
  "grain_boundary_misorientation_distribution.csv";
  "grain_boundary_misorientation_distribution.png";
  "grain_boundary_misorientation_metrics.png";
  "grain_boundary_misorientation_summary.csv"
];
for artifactName = artifactNames'
  artifactPath = fullfile(outputDir, artifactName);
  if isfile(artifactPath)
    delete(artifactPath);
  end
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
distributionEdgesDeg = 1:0.5:94;
networkEstimator = ...
  "MTEX thresholded Ti-Hex/Ti-Hex boundary-segment network";
eligibleSources = "innerBoundary LAGB + grains.boundary HAGB";

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
  assert(~contains(inputFile, "_denoised.ctf"), ...
    "Denoised CTF files are not valid inputs: %s", inputFile);
  fprintf("IMPORT=%s\n", inputFile);

  ebsdFull = EBSD.load(inputFile, 'convertEuler2SpatialReferenceFrame');
  scanAudit = audit_full_native_scan(ebsdFull, sampleNames(sampleIndex));
  mappedEbsdId = double(ebsdFull.id(:));
  mappedXY = [double(ebsdFull.x(:)), double(ebsdFull.y(:))];

  for floorIndex = 1:floorCount
    floorDeg = detectionFloorsDeg(floorIndex);
    [grains, ebsdFull.grainId] = calcGrains(ebsdFull, 'unitCell', ...
      'threshold', [floorDeg classificationAngleDeg] * degree, ...
      'minPixel', minPixel);

    inner = grains.innerBoundary("Ti-Hex", "Ti-Hex");
    outer = grains.boundary("Ti-Hex", "Ti-Hex");
    innerThetaDeg = double(angle(inner.misorientation) / degree);
    outerThetaDeg = double(angle(outer.misorientation) / degree);
    innerLengthUm = double(inner.segLength(:));
    outerLengthUm = double(outer.segLength(:));

    [thetaDeg, segLengthUm, sourceAudit] = ...
      partition_ti_hex_boundary_segments(innerThetaDeg, innerLengthUm, ...
      outerThetaDeg, outerLengthUm, floorDeg, classificationAngleDeg);
    assert(~isempty(thetaDeg), ...
      "No eligible Ti-Hex boundaries for %s at %.1f degrees.", ...
      sampleNames(sampleIndex), floorDeg);

    sourceEndpointIds = double([inner.ebsdId; outer.ebsdId]);
    topologyAudit = audit_native_grid_pairs(sourceEndpointIds, ...
      mappedEbsdId, mappedXY, scanAudit.native_x_step_um, ...
      scanAudit.native_y_step_um, ...
      sampleNames(sampleIndex) + " Ti-Hex/Ti-Hex source faces");
    assert(topologyAudit.endpoint_pair_count == ...
      sourceAudit.source_network_segment_count, ...
      "Source segments and native endpoint pairs do not match.");

    boundaryStats = calculate_boundary_metrics(thetaDeg, segLengthUm, ...
      classificationAngleDeg, scanAudit.indexed_ti_hex_area_um2);
    assert(boundaryStats.eligible_segment_count == ...
      sourceAudit.eligible_segment_count);
    assert(abs(boundaryStats.total_eligible_boundary_length_um - ...
      sourceAudit.total_eligible_boundary_length_um) < 1e-8);

    row = struct();
    row.sample = sampleNames(sampleIndex);
    row.folder = folderNames(sampleIndex);
    row.input_file = inputNames(sampleIndex);
    row.cold_reduction_percent = coldReductionPercent(sampleIndex);
    row.detection_floor_deg = floorDeg;
    row.classification_angle_deg = classificationAngleDeg;
    row.min_pixel = minPixel;
    row.boundary_network_estimator = networkEstimator;
    row.eligible_boundary_sources = eligibleSources;
    row.scan_unit = scanAudit.scan_unit;
    row.native_x_cell_count = scanAudit.native_x_cell_count;
    row.native_y_cell_count = scanAudit.native_y_cell_count;
    row.native_x_step_um = scanAudit.native_x_step_um;
    row.native_y_step_um = scanAudit.native_y_step_um;
    row.pixel_area_um2 = scanAudit.pixel_area_um2;
    row.mapped_pixel_count = scanAudit.mapped_pixel_count;
    row.mapped_area_um2 = scanAudit.mapped_area_um2;
    row.unindexed_pixel_count = scanAudit.unindexed_pixel_count;
    row.indexed_ti_hex_pixel_count = ...
      scanAudit.indexed_ti_hex_pixel_count;
    row.indexed_ti_hex_area_um2 = scanAudit.indexed_ti_hex_area_um2;
    row.indexed_ti_hex_fraction_of_mapped_area = ...
      scanAudit.indexed_ti_hex_fraction_of_mapped_area;
    row.ti_hex_grain_count = length(grains("Ti-Hex"));
    row.inner_ti_ti_segment_count = ...
      sourceAudit.inner_ti_ti_segment_count;
    row.inner_ti_ti_boundary_length_um = ...
      sourceAudit.inner_ti_ti_boundary_length_um;
    row.outer_ti_ti_segment_count = ...
      sourceAudit.outer_ti_ti_segment_count;
    row.outer_ti_ti_boundary_length_um = ...
      sourceAudit.outer_ti_ti_boundary_length_um;
    row.audit_inner_lagb_segment_count = ...
      sourceAudit.audit_inner_lagb_segment_count;
    row.audit_inner_lagb_length_um = ...
      sourceAudit.audit_inner_lagb_length_um;
    row.audit_inner_high_angle_segment_count = ...
      sourceAudit.audit_inner_high_angle_segment_count;
    row.audit_inner_high_angle_length_um = ...
      sourceAudit.audit_inner_high_angle_length_um;
    row.audit_outer_low_angle_segment_count = ...
      sourceAudit.audit_outer_low_angle_segment_count;
    row.audit_outer_low_angle_length_um = ...
      sourceAudit.audit_outer_low_angle_length_um;
    row.audit_outer_hagb_segment_count = ...
      sourceAudit.audit_outer_hagb_segment_count;
    row.audit_outer_hagb_length_um = ...
      sourceAudit.audit_outer_hagb_length_um;
    row.source_network_segment_count = ...
      sourceAudit.source_network_segment_count;
    row.source_network_boundary_length_um = ...
      sourceAudit.source_network_boundary_length_um;
    row.excluded_cross_class_segment_count = ...
      sourceAudit.excluded_cross_class_segment_count;
    row.excluded_cross_class_boundary_length_um = ...
      sourceAudit.excluded_cross_class_boundary_length_um;
    row.native_ti_ti_endpoint_pair_count = ...
      topologyAudit.endpoint_pair_count;
    row.nonlocal_ti_ti_endpoint_pair_count = ...
      topologyAudit.nonlocal_endpoint_pair_count;
    row.max_ti_ti_endpoint_distance_um = ...
      topologyAudit.max_endpoint_distance_um;
    row.eligible_segment_count = sourceAudit.eligible_segment_count;
    row.total_eligible_boundary_length_um = ...
      sourceAudit.total_eligible_boundary_length_um;
    row.lagb_segment_count = boundaryStats.lagb_segment_count;
    row.hagb_segment_count = boundaryStats.hagb_segment_count;
    row.lagb_length_um = boundaryStats.lagb_length_um;
    row.hagb_length_um = boundaryStats.hagb_length_um;
    row.lagb_length_fraction = boundaryStats.lagb_length_fraction;
    row.hagb_length_fraction = boundaryStats.hagb_length_fraction;
    row.total_eligible_length_density_um_per_um2 = ...
      boundaryStats.total_eligible_length_density_um_per_um2;
    row.lagb_length_density_um_per_um2 = ...
      boundaryStats.lagb_length_density_um_per_um2;
    row.hagb_length_density_um_per_um2 = ...
      boundaryStats.hagb_length_density_um_per_um2;
    row.weighted_mean_angle_deg = boundaryStats.weighted_mean_angle_deg;
    row.weighted_median_angle_deg = ...
      boundaryStats.weighted_median_angle_deg;
    row.min_angle_deg = boundaryStats.min_angle_deg;
    row.max_angle_deg = boundaryStats.max_angle_deg;

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
          intervalLength / sourceAudit.total_eligible_boundary_length_um;
      end
      summaryRows{sampleIndex} = struct2table(summaryRow);

      [weightedStats, sampleDistribution] = ...
        summarize_weighted_boundary_angles(thetaDeg, segLengthUm, ...
        distributionEdgesDeg);
      assert(abs(weightedStats.total_boundary_length_um - ...
        sourceAudit.total_eligible_boundary_length_um) < 1e-8);
      blockHeight = height(sampleDistribution);
      distributionBlocks{sampleIndex} = addvars(sampleDistribution, ...
        repmat(sampleNames(sampleIndex), blockHeight, 1), ...
        repmat(folderNames(sampleIndex), blockHeight, 1), ...
        repmat(inputNames(sampleIndex), blockHeight, 1), ...
        repmat(coldReductionPercent(sampleIndex), blockHeight, 1), ...
        repmat(floorDeg, blockHeight, 1), ...
        repmat(classificationAngleDeg, blockHeight, 1), ...
        repmat(minPixel, blockHeight, 1), ...
        repmat(networkEstimator, blockHeight, 1), ...
        repmat(eligibleSources, blockHeight, 1), ...
        repmat(scanAudit.scan_unit, blockHeight, 1), ...
        repmat(scanAudit.indexed_ti_hex_area_um2, blockHeight, 1), ...
        'Before', 1, 'NewVariableNames', ...
        {'sample', 'folder', 'input_file', 'cold_reduction_percent', ...
        'detection_floor_deg', 'classification_angle_deg', 'min_pixel', ...
        'boundary_network_estimator', 'eligible_boundary_sources', ...
        'scan_unit', 'indexed_ti_hex_area_um2'});
    end

    fprintf("SAMPLE=%s FLOOR_DEG=%.1f TI_GRAINS=%d SOURCE_UM=%.6f " + ...
      "ELIGIBLE_UM=%.6f INNER_HIGH_UM=%.6f OUTER_LOW_UM=%.6f " + ...
      "LAGB_FRACTION=%.6f\n", sampleNames(sampleIndex), floorDeg, ...
      row.ti_hex_grain_count, row.source_network_boundary_length_um, ...
      row.total_eligible_boundary_length_um, ...
      row.audit_inner_high_angle_length_um, ...
      row.audit_outer_low_angle_length_um, row.lagb_length_fraction);
    clear grains inner outer innerThetaDeg outerThetaDeg
    clear innerLengthUm outerLengthUm thetaDeg segLengthUm
  end
  clear ebsdFull
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
plot_boundary_metrics(sensitivity, metricsFigureFile);

fprintf("SUMMARY=%s\n", summaryFile);
fprintf("SENSITIVITY=%s\n", sensitivityFile);
fprintf("DISTRIBUTION=%s\n", distributionFile);
fprintf("DISTRIBUTION_FIGURE=%s\n", distributionFigureFile);
fprintf("METRICS_FIGURE=%s\n", metricsFigureFile);
end

function audit = audit_full_native_scan(ebsdFull, sampleName)
coordinateTolerance = 1e-10;
xValues = unique(double(ebsdFull.x(:)));
yValues = unique(double(ebsdFull.y(:)));
assert(numel(xValues) > 1 && numel(yValues) > 1, ...
  "Native scan dimensions are invalid for %s.", sampleName);
xSteps = diff(xValues);
ySteps = diff(yValues);
xStepUm = median(xSteps);
yStepUm = median(ySteps);
pixelAreaUm2 = double(polyarea(ebsdFull.unitCell.x, ebsdFull.unitCell.y));
mappedPixelCount = length(ebsdFull);
unindexedPixelCount = length(ebsdFull("notIndexed"));
indexedTiHexPixelCount = length(ebsdFull("Ti-Hex"));
mappedAreaUm2 = mappedPixelCount * pixelAreaUm2;
indexedTiHexAreaUm2 = indexedTiHexPixelCount * pixelAreaUm2;
mappedEbsdId = double(ebsdFull.id(:));

assert(string(ebsdFull.scanUnit) == "um", ...
  "Expected micrometre scan units for %s.", sampleName);
assert(numel(xValues) == 600 && numel(yValues) == 600, ...
  "Expected a native 600-by-600 grid for %s.", sampleName);
assert(mappedPixelCount == 360000 && ...
  mappedPixelCount == numel(xValues) * numel(yValues), ...
  "Expected all 360000 native mapped sites for %s.", sampleName);
assert(all(abs(xSteps - 0.5) < coordinateTolerance) && ...
  all(abs(ySteps - 0.5) < coordinateTolerance), ...
  "Expected a 0.5 um native step for %s.", sampleName);
assert(abs(pixelAreaUm2 - 0.25) < coordinateTolerance && ...
  abs(mappedAreaUm2 - 90000) < coordinateTolerance, ...
  "Native pixel or mapped area is invalid for %s.", sampleName);
assert(unindexedPixelCount > 0 && ...
  indexedTiHexPixelCount <= mappedPixelCount - unindexedPixelCount, ...
  "Indexed and unindexed counts are inconsistent for %s.", sampleName);
assert(numel(unique(mappedEbsdId)) == mappedPixelCount, ...
  "Persistent EBSD IDs are not unique for %s.", sampleName);

audit.scan_unit = string(ebsdFull.scanUnit);
audit.native_x_cell_count = numel(xValues);
audit.native_y_cell_count = numel(yValues);
audit.native_x_step_um = xStepUm;
audit.native_y_step_um = yStepUm;
audit.pixel_area_um2 = pixelAreaUm2;
audit.mapped_pixel_count = mappedPixelCount;
audit.mapped_area_um2 = mappedAreaUm2;
audit.unindexed_pixel_count = unindexedPixelCount;
audit.indexed_ti_hex_pixel_count = indexedTiHexPixelCount;
audit.indexed_ti_hex_area_um2 = indexedTiHexAreaUm2;
audit.indexed_ti_hex_fraction_of_mapped_area = ...
  indexedTiHexAreaUm2 / mappedAreaUm2;
end

function stats = calculate_boundary_metrics(thetaDeg, segLengthUm, ...
  classificationAngleDeg, indexedAreaUm2)
assert(~isempty(thetaDeg) && numel(thetaDeg) == numel(segLengthUm));
assert(all(isfinite(thetaDeg)) && ...
  all(isfinite(segLengthUm) & segLengthUm > 0));
assert(isfinite(indexedAreaUm2) && indexedAreaUm2 > 0);
totalLengthUm = sum(segLengthUm);
lagbMask = thetaDeg < classificationAngleDeg;
hagbMask = thetaDeg >= classificationAngleDeg;
assert(all(xor(lagbMask, hagbMask)));
lagbLengthUm = sum(segLengthUm(lagbMask));
hagbLengthUm = sum(segLengthUm(hagbMask));
[sortedTheta, order] = sort(thetaDeg);
sortedLength = segLengthUm(order);
weightedMedianIndex = find(cumsum(sortedLength) >= 0.5 * totalLengthUm, ...
  1);

stats.eligible_segment_count = numel(thetaDeg);
stats.total_eligible_boundary_length_um = totalLengthUm;
stats.lagb_segment_count = nnz(lagbMask);
stats.hagb_segment_count = nnz(hagbMask);
stats.lagb_length_um = lagbLengthUm;
stats.hagb_length_um = hagbLengthUm;
stats.lagb_length_fraction = lagbLengthUm / totalLengthUm;
stats.hagb_length_fraction = hagbLengthUm / totalLengthUm;
stats.total_eligible_length_density_um_per_um2 = ...
  totalLengthUm / indexedAreaUm2;
stats.lagb_length_density_um_per_um2 = lagbLengthUm / indexedAreaUm2;
stats.hagb_length_density_um_per_um2 = hagbLengthUm / indexedAreaUm2;
stats.weighted_mean_angle_deg = ...
  sum(thetaDeg .* segLengthUm) / totalLengthUm;
stats.weighted_median_angle_deg = sortedTheta(weightedMedianIndex);
stats.min_angle_deg = min(thetaDeg);
stats.max_angle_deg = max(thetaDeg);
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
    "Ti-Hex boundary-segment misorientation (degree)");
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
title(layout, ...
  "Primary detection floor = 1 degree (0.5-degree bins; no smoothing)");

export_clean_png(figureHandle, outputFile);
clear cleanupFigure
end

function plot_boundary_metrics(sensitivity, outputFile)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Position", [100, 100, 1200, 800]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 2, 1, "Padding", "compact", ...
  "TileSpacing", "compact");
detectionFloorsDeg = unique(sensitivity.detection_floor_deg, "stable");
colors = lines(numel(detectionFloorsDeg));
markers = ["o"; "s"; "^"];
metricNames = ["lagb_length_fraction"; ...
  "lagb_length_density_um_per_um2"];
yLabels = ["LAGB length fraction (%)"; ...
  "LAGB length density (\mum \mum^{-2})"];
panelTitles = ["Low-angle boundary-segment length fraction"; ...
  "Low-angle boundary-segment length density"];

for panelIndex = 1:2
  axesHandle = nexttile(layout);
  hold(axesHandle, "on");
  for floorIndex = 1:numel(detectionFloorsDeg)
    floorDeg = detectionFloorsDeg(floorIndex);
    rows = sensitivity.detection_floor_deg == floorDeg;
    xValues = sensitivity.cold_reduction_percent(rows);
    yValues = sensitivity.(metricNames(panelIndex))(rows);
    if panelIndex == 1
      yValues = 100 * yValues;
    end
    plot(axesHandle, xValues, yValues, ...
      "Color", colors(floorIndex,:), "LineWidth", 1.4, ...
      "Marker", markers(floorIndex), "MarkerSize", 6, ...
      "MarkerFaceColor", colors(floorIndex,:), ...
      "DisplayName", sprintf("Detection floor = %.1f degree%s", ...
      floorDeg, primary_suffix(floorDeg)));
  end

  primaryRows = sensitivity.detection_floor_deg == 1;
  primaryX = sensitivity.cold_reduction_percent(primaryRows);
  primaryY = sensitivity.(metricNames(panelIndex))(primaryRows);
  if panelIndex == 1
    primaryY = 100 * primaryY;
  end
  for sampleIndex = 1:nnz(primaryRows)
    primarySampleNames = sensitivity.sample(primaryRows);
    text(axesHandle, primaryX(sampleIndex), primaryY(sampleIndex), ...
      "  " + primarySampleNames(sampleIndex), ...
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
  if panelIndex == 1
    legend(axesHandle, "Location", "eastoutside");
  end
end
title(layout, "Detection-floor sensitivity (no smoothing)");

export_clean_png(figureHandle, outputFile);
clear cleanupFigure
end

function suffix = primary_suffix(floorDeg)
if floorDeg == 1
  suffix = " (primary)";
else
  suffix = "";
end
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

function validate_tables(summary, sensitivity, distribution, sampleNames, ...
  detectionFloorsDeg, distributionEdgesDeg)
countTolerance = 0;
lengthTolerance = 1e-8;
fractionTolerance = 1e-10;
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

assert(all(sensitivity.scan_unit == "um"));
assert(all(sensitivity.native_x_cell_count == 600 & ...
  sensitivity.native_y_cell_count == 600));
assert(all(abs(sensitivity.native_x_step_um - 0.5) < fractionTolerance & ...
  abs(sensitivity.native_y_step_um - 0.5) < fractionTolerance));
assert(all(sensitivity.mapped_pixel_count == 360000));
assert(all(abs(sensitivity.pixel_area_um2 - 0.25) < fractionTolerance));
assert(all(abs(sensitivity.mapped_area_um2 - 90000) < lengthTolerance));
assert(all(sensitivity.unindexed_pixel_count > 0));
assert(all(sensitivity.indexed_ti_hex_pixel_count <= ...
  sensitivity.mapped_pixel_count - sensitivity.unindexed_pixel_count));
assert(all(abs(sensitivity.indexed_ti_hex_area_um2 - ...
  sensitivity.indexed_ti_hex_pixel_count .* ...
  sensitivity.pixel_area_um2) < lengthTolerance));
assert(all(abs(sensitivity.indexed_ti_hex_fraction_of_mapped_area - ...
  sensitivity.indexed_ti_hex_area_um2 ./ sensitivity.mapped_area_um2) ...
  < fractionTolerance));
assert(all(sensitivity.ti_hex_grain_count > 0 & ...
  sensitivity.ti_hex_grain_count == fix(sensitivity.ti_hex_grain_count)));

assert(all(sensitivity.inner_ti_ti_segment_count == ...
  sensitivity.audit_inner_lagb_segment_count + ...
  sensitivity.audit_inner_high_angle_segment_count));
assert(all(sensitivity.outer_ti_ti_segment_count == ...
  sensitivity.audit_outer_low_angle_segment_count + ...
  sensitivity.audit_outer_hagb_segment_count));
assert(all(sensitivity.source_network_segment_count == ...
  sensitivity.inner_ti_ti_segment_count + ...
  sensitivity.outer_ti_ti_segment_count));
assert(all(sensitivity.eligible_segment_count == ...
  sensitivity.audit_inner_lagb_segment_count + ...
  sensitivity.audit_outer_hagb_segment_count));
assert(all(sensitivity.excluded_cross_class_segment_count == ...
  sensitivity.audit_inner_high_angle_segment_count + ...
  sensitivity.audit_outer_low_angle_segment_count));
assert(all(sensitivity.source_network_segment_count == ...
  sensitivity.eligible_segment_count + ...
  sensitivity.excluded_cross_class_segment_count));

assert(all(abs(sensitivity.inner_ti_ti_boundary_length_um - ...
  sensitivity.audit_inner_lagb_length_um - ...
  sensitivity.audit_inner_high_angle_length_um) < lengthTolerance));
assert(all(abs(sensitivity.outer_ti_ti_boundary_length_um - ...
  sensitivity.audit_outer_low_angle_length_um - ...
  sensitivity.audit_outer_hagb_length_um) < lengthTolerance));
assert(all(abs(sensitivity.source_network_boundary_length_um - ...
  sensitivity.inner_ti_ti_boundary_length_um - ...
  sensitivity.outer_ti_ti_boundary_length_um) < lengthTolerance));
assert(all(abs(sensitivity.total_eligible_boundary_length_um - ...
  sensitivity.audit_inner_lagb_length_um - ...
  sensitivity.audit_outer_hagb_length_um) < lengthTolerance));
assert(all(abs(sensitivity.excluded_cross_class_boundary_length_um - ...
  sensitivity.audit_inner_high_angle_length_um - ...
  sensitivity.audit_outer_low_angle_length_um) < lengthTolerance));
assert(all(abs(sensitivity.source_network_boundary_length_um - ...
  sensitivity.total_eligible_boundary_length_um - ...
  sensitivity.excluded_cross_class_boundary_length_um) < lengthTolerance));

assert(all(sensitivity.native_ti_ti_endpoint_pair_count == ...
  sensitivity.source_network_segment_count));
assert(all(sensitivity.nonlocal_ti_ti_endpoint_pair_count == countTolerance));
assert(all(sensitivity.max_ti_ti_endpoint_distance_um <= 0.5 + ...
  fractionTolerance));
assert(all(sensitivity.lagb_segment_count == ...
  sensitivity.audit_inner_lagb_segment_count));
assert(all(sensitivity.hagb_segment_count == ...
  sensitivity.audit_outer_hagb_segment_count));
assert(all(abs(sensitivity.lagb_length_um - ...
  sensitivity.audit_inner_lagb_length_um) < lengthTolerance));
assert(all(abs(sensitivity.hagb_length_um - ...
  sensitivity.audit_outer_hagb_length_um) < lengthTolerance));
assert(all(abs(sensitivity.lagb_length_um + sensitivity.hagb_length_um - ...
  sensitivity.total_eligible_boundary_length_um) < lengthTolerance));
assert(all(abs(sensitivity.lagb_length_fraction + ...
  sensitivity.hagb_length_fraction - 1) < fractionTolerance));
assert(all(abs(sensitivity.total_eligible_length_density_um_per_um2 - ...
  sensitivity.total_eligible_boundary_length_um ./ ...
  sensitivity.indexed_ti_hex_area_um2) < fractionTolerance));

primaryRows = sensitivity.detection_floor_deg == 1;
baseNames = sensitivity.Properties.VariableNames;
assert(isequal(summary(:, baseNames), sensitivity(primaryRows, :)), ...
  "Primary summary does not reproduce the 1-degree sensitivity rows.");
intervalLength = summary.length_1_2_um + summary.length_2_5_um + ...
  summary.length_5_10_um + summary.length_10_15_um + ...
  summary.length_15_94_um;
assert(all(abs(intervalLength - ...
  summary.total_eligible_boundary_length_um) < lengthTolerance));
intervalFraction = summary.fraction_1_2 + summary.fraction_2_5 + ...
  summary.fraction_5_10 + summary.fraction_10_15 + ...
  summary.fraction_15_94;
assert(all(abs(intervalFraction - 1) < fractionTolerance));
assert(all(abs(summary.length_15_94_um - summary.hagb_length_um) ...
  < lengthTolerance));

for sampleIndex = 1:numel(sampleNames)
  sampleName = sampleNames(sampleIndex);
  sensitivityMask = sensitivity.sample == sampleName;
  assert(nnz(sensitivityMask) == numel(detectionFloorsDeg));
  assert(isequal(sort(sensitivity.detection_floor_deg(sensitivityMask)), ...
    detectionFloorsDeg));
  summaryMask = summary.sample == sampleName;
  distributionMask = distribution.sample == sampleName;
  assert(nnz(summaryMask) == 1);
  assert(nnz(distributionMask) == numel(distributionEdgesDeg) - 1);
  assert(all(distribution.detection_floor_deg(distributionMask) == 1));
  assert(sum(distribution.segment_count(distributionMask)) == ...
    summary.eligible_segment_count(summaryMask));
  assert(abs(sum(distribution.boundary_length_um(distributionMask)) - ...
    summary.total_eligible_boundary_length_um(summaryMask)) ...
    < lengthTolerance);
  assert(abs(sum(distribution.length_fraction(distributionMask)) - 1) ...
    < fractionTolerance);
  assert(abs(sum(distribution.probability_density_per_degree( ...
    distributionMask) .* distribution.bin_width_deg(distributionMask)) ...
    - 1) < fractionTolerance);
  distributionIndices = find(distributionMask);
  assert(abs(distribution.cumulative_length_fraction( ...
    distributionIndices(end)) - 1) < fractionTolerance);
end
end
