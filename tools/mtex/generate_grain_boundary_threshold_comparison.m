function [summary, distribution] = ...
  generate_grain_boundary_threshold_comparison(scanRoot, outputDir)
%GENERATE_GRAIN_BOUNDARY_THRESHOLD_COMPARISON Compare 5 and 15 degree cutoffs.

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
  "grain_boundary_threshold_comparison_summary.csv";
  "grain_boundary_threshold_comparison_distribution.csv"
];
for artifactName = artifactNames'
  artifactPath = fullfile(outputDir, artifactName);
  if isfile(artifactPath)
    delete(artifactPath);
  end
end

folderNames = ["d7";"d6_48";"d6_02";"d5_6";"d5_25";"d5"];
inputNames = [
  "ebsd_sample_7_map_15.ctf";
  "ebsd_sample_648_map_13.ctf";
  "ebsd_sample_602_map_11.ctf";
  "ebsd_sample_56_map_9.ctf";
  "ebsd_sample_525_map_7.ctf";
  "ebsd_sample_5_map_3.ctf"
];
sampleNames = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"];
coldReductionPercent = [0;14.31;26.04;36.00;43.75;48.98];
detectionFloorDeg = 2;
reconstructionUpperAngleDeg = 15;
classificationAnglesDeg = [5;15];
minPixel = 5;
distributionEdgesDeg = 1:0.5:94;
networkEstimator = ...
  "MTEX thresholded Ti-Hex/Ti-Hex boundary-segment network";
eligibleSources = "innerBoundary LAGB + grains.boundary HAGB";

sampleCount = numel(sampleNames);
thresholdCount = numel(classificationAnglesDeg);
summaryRows = cell(sampleCount * thresholdCount, 1);
distributionBlocks = cell(sampleCount * thresholdCount, 1);
rowIndex = 0;

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

  [grains, ebsdFull.grainId] = calcGrains(ebsdFull, 'unitCell', ...
    'threshold', [2 15] * degree, 'minPixel', minPixel);
  inner = grains.innerBoundary("Ti-Hex", "Ti-Hex");
  outer = grains.boundary("Ti-Hex", "Ti-Hex");
  innerThetaDeg = double(angle(inner.misorientation) / degree);
  outerThetaDeg = double(angle(outer.misorientation) / degree);
  innerLengthUm = double(inner.segLength(:));
  outerLengthUm = double(outer.segLength(:));

  [thetaDeg, segLengthUm, sourceAudit] = ...
    partition_ti_hex_boundary_segments(innerThetaDeg, innerLengthUm, ...
    outerThetaDeg, outerLengthUm, detectionFloorDeg, ...
    reconstructionUpperAngleDeg);
  assert(~isempty(thetaDeg), "No eligible Ti-Hex boundaries for %s.", ...
    sampleNames(sampleIndex));

  sourceEndpointIds = double([inner.ebsdId;outer.ebsdId]);
  topologyAudit = audit_native_grid_pairs(sourceEndpointIds, mappedEbsdId, ...
    mappedXY, scanAudit.native_x_step_um, scanAudit.native_y_step_um, ...
    sampleNames(sampleIndex) + " Ti-Hex/Ti-Hex source faces");
  assert(topologyAudit.endpoint_pair_count == ...
    sourceAudit.source_network_segment_count);

  [weightedStats, sampleDistribution] = ...
    summarize_weighted_boundary_angles(thetaDeg, segLengthUm, ...
    distributionEdgesDeg);
  assert(abs(weightedStats.total_boundary_length_um - ...
    sourceAudit.total_eligible_boundary_length_um) < 1e-8);

  for thresholdIndex = 1:thresholdCount
    classificationAngleDeg = classificationAnglesDeg(thresholdIndex);
    metrics = calculate_boundary_threshold_metrics(thetaDeg, segLengthUm, ...
      detectionFloorDeg, classificationAngleDeg, ...
      scanAudit.indexed_ti_hex_area_um2);
    rowIndex = rowIndex + 1;

    row = struct();
    row.sample = sampleNames(sampleIndex);
    row.folder = folderNames(sampleIndex);
    row.input_file = inputNames(sampleIndex);
    row.cold_reduction_percent = coldReductionPercent(sampleIndex);
    row.detection_floor_deg = detectionFloorDeg;
    row.reconstruction_upper_angle_deg = reconstructionUpperAngleDeg;
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
    row.indexed_ti_hex_pixel_count = scanAudit.indexed_ti_hex_pixel_count;
    row.indexed_ti_hex_area_um2 = scanAudit.indexed_ti_hex_area_um2;
    row.ti_hex_grain_count = length(grains("Ti-Hex"));
    row.source_network_segment_count = ...
      sourceAudit.source_network_segment_count;
    row.source_network_boundary_length_um = ...
      sourceAudit.source_network_boundary_length_um;
    row.eligible_segment_count = metrics.eligible_segment_count;
    row.total_eligible_boundary_length_um = ...
      metrics.total_eligible_boundary_length_um;
    row.native_ti_ti_endpoint_pair_count = ...
      topologyAudit.endpoint_pair_count;
    row.nonlocal_ti_ti_endpoint_pair_count = ...
      topologyAudit.nonlocal_endpoint_pair_count;
    row.max_ti_ti_endpoint_distance_um = ...
      topologyAudit.max_endpoint_distance_um;
    row.lagb_segment_count = metrics.lagb_segment_count;
    row.hagb_segment_count = metrics.hagb_segment_count;
    row.lagb_length_um = metrics.lagb_length_um;
    row.hagb_length_um = metrics.hagb_length_um;
    row.lagb_length_fraction = metrics.lagb_length_fraction;
    row.hagb_length_fraction = metrics.hagb_length_fraction;
    row.total_eligible_length_density_um_per_um2 = ...
      metrics.total_eligible_length_density_um_per_um2;
    row.lagb_length_density_um_per_um2 = ...
      metrics.lagb_length_density_um_per_um2;
    row.hagb_length_density_um_per_um2 = ...
      metrics.hagb_length_density_um_per_um2;
    row.weighted_mean_angle_deg = metrics.weighted_mean_angle_deg;
    row.weighted_median_angle_deg = metrics.weighted_median_angle_deg;
    row.min_angle_deg = metrics.min_angle_deg;
    row.max_angle_deg = metrics.max_angle_deg;
    summaryRows{rowIndex} = struct2table(row);

    blockHeight = height(sampleDistribution);
    distributionBlocks{rowIndex} = addvars(sampleDistribution, ...
      repmat(sampleNames(sampleIndex), blockHeight, 1), ...
      repmat(folderNames(sampleIndex), blockHeight, 1), ...
      repmat(inputNames(sampleIndex), blockHeight, 1), ...
      repmat(coldReductionPercent(sampleIndex), blockHeight, 1), ...
      repmat(detectionFloorDeg, blockHeight, 1), ...
      repmat(reconstructionUpperAngleDeg, blockHeight, 1), ...
      repmat(classificationAngleDeg, blockHeight, 1), ...
      'Before', 1, 'NewVariableNames', ...
      {'sample','folder','input_file','cold_reduction_percent', ...
      'detection_floor_deg','reconstruction_upper_angle_deg', ...
      'classification_angle_deg'});

    fprintf("SAMPLE=%s CUTOFF_DEG=%g LAGB_FRACTION=%.6f " + ...
      "LAGB_DENSITY=%.6f\n", sampleNames(sampleIndex), ...
      classificationAngleDeg, metrics.lagb_length_fraction, ...
      metrics.lagb_length_density_um_per_um2);
  end
  clear ebsdFull grains inner outer
end

summary = vertcat(summaryRows{:});
distribution = vertcat(distributionBlocks{:});
validate_comparison_tables(summary, distribution, sampleNames, ...
  classificationAnglesDeg);

summaryFile = fullfile(outputDir, ...
  "grain_boundary_threshold_comparison_summary.csv");
distributionFile = fullfile(outputDir, ...
  "grain_boundary_threshold_comparison_distribution.csv");
writetable(summary, summaryFile);
writetable(distribution, distributionFile);
fprintf("SUMMARY=%s\n", summaryFile);
fprintf("DISTRIBUTION=%s\n", distributionFile);
end

function audit = audit_full_native_scan(ebsdFull, ~)
coordinateTolerance = 1e-10;
xValues = unique(double(ebsdFull.x(:)));
yValues = unique(double(ebsdFull.y(:)));
assert(numel(xValues) > 1 && numel(yValues) > 1);
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

assert(string(ebsdFull.scanUnit) == "um");
assert(numel(xValues) == 600 && numel(yValues) == 600);
assert(mappedPixelCount == 360000);
assert(all(abs(xSteps - 0.5) < coordinateTolerance));
assert(all(abs(ySteps - 0.5) < coordinateTolerance));
assert(abs(pixelAreaUm2 - 0.25) < coordinateTolerance);
assert(abs(mappedAreaUm2 - 90000) < coordinateTolerance);
assert(unindexedPixelCount > 0);

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
end

function validate_comparison_tables(summary, distribution, sampleNames, ...
  classificationAnglesDeg)
assert(height(summary) == numel(sampleNames) * ...
  numel(classificationAnglesDeg));
assert(height(distribution) == height(summary) * 186);
assert(all(summary.detection_floor_deg == 2));
assert(all(summary.reconstruction_upper_angle_deg == 15));
assert(all(summary.nonlocal_ti_ti_endpoint_pair_count == 0));
assert(all(abs(summary.lagb_length_um + summary.hagb_length_um - ...
  summary.total_eligible_boundary_length_um) < 1e-8));
assert(all(abs(summary.lagb_length_fraction + ...
  summary.hagb_length_fraction - 1) < 1e-10));
for sampleName = sampleNames'
  rows = summary.sample == sampleName;
  assert(nnz(rows) == numel(classificationAnglesDeg));
  assert(isequal(summary.classification_angle_deg(rows), ...
    classificationAnglesDeg));
  assert(max(summary.total_eligible_boundary_length_um(rows)) - ...
    min(summary.total_eligible_boundary_length_um(rows)) < 1e-8);
end
end
