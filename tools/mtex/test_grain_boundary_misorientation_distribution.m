function test_grain_boundary_misorientation_distribution(scanRoot, outputDir)
arguments
  scanRoot (1,1) string = ""
  outputDir (1,1) string = ""
end

innerThetaDeg = [1.2;15;20];
innerLengthUm = [2;3;4];
outerThetaDeg = [1.1;14.9;15;30];
outerLengthUm = [5;6;7;8];
[thetaDeg, segLengthUm, sourceAudit, sourceMasks] = ...
  partition_ti_hex_boundary_segments(innerThetaDeg, innerLengthUm, ...
  outerThetaDeg, outerLengthUm, 1, 15);
assert(isequal(thetaDeg, [1.2;15;30]));
assert(isequal(segLengthUm, [2;7;8]));
assert(isequal(sourceMasks.inner_lagb, logical([1;0;0])));
assert(isequal(sourceMasks.inner_high_angle, logical([0;1;1])));
assert(isequal(sourceMasks.outer_low_angle, logical([1;1;0;0])));
assert(isequal(sourceMasks.outer_hagb, logical([0;0;1;1])));
assert(sourceAudit.audit_inner_lagb_segment_count == 1);
assert(sourceAudit.audit_inner_lagb_length_um == 2);
assert(sourceAudit.audit_inner_high_angle_segment_count == 2);
assert(sourceAudit.audit_inner_high_angle_length_um == 7);
assert(sourceAudit.audit_outer_low_angle_segment_count == 2);
assert(sourceAudit.audit_outer_low_angle_length_um == 11);
assert(sourceAudit.audit_outer_hagb_segment_count == 2);
assert(sourceAudit.audit_outer_hagb_length_um == 15);
assert(sourceAudit.source_network_segment_count == 7);
assert(sourceAudit.source_network_boundary_length_um == 35);
assert(sourceAudit.eligible_segment_count == 3);
assert(sourceAudit.total_eligible_boundary_length_um == 17);
assert(sourceAudit.excluded_cross_class_segment_count == 4);
assert(sourceAudit.excluded_cross_class_boundary_length_um == 18);

mappedEbsdId = [305;101;999;707];
mappedXY = [0.5,0;0,0;0.5,0.5;0,0.5];
pairAudit = audit_native_grid_pairs([101,305;305,999], ...
  mappedEbsdId, mappedXY, 0.5, 0.5, "synthetic local pairs");
assert(pairAudit.endpoint_pair_count == 2);
assert(pairAudit.nonlocal_endpoint_pair_count == 0);
assert(abs(pairAudit.max_endpoint_distance_um - 0.5) < 1e-12);
didRejectNonlocalPair = false;
try
  audit_native_grid_pairs([101,999], mappedEbsdId, mappedXY, 0.5, 0.5, ...
    "synthetic nonlocal pair");
catch exception
  didRejectNonlocalPair = ...
    exception.identifier == "audit_native_grid_pairs:NonlocalEndpoint";
end
assert(didRejectNonlocalPair, ...
  "The native-grid audit did not reject a non-face-neighbor pair.");
didRejectUnmappedId = false;
try
  audit_native_grid_pairs([101,777], mappedEbsdId, mappedXY, 0.5, 0.5, ...
    "synthetic unmapped endpoint");
catch exception
  didRejectUnmappedId = ...
    exception.identifier == "audit_native_grid_pairs:UnmappedEndpoint";
end
assert(didRejectUnmappedId, ...
  "The native-grid audit did not reject an unmapped persistent EBSD ID.");
didRejectIncompleteGrid = false;
try
  duplicateXY = [0.5,0;0,0;0.5,0.5;0.5,0.5];
  audit_native_grid_pairs([101,305], mappedEbsdId, duplicateXY, ...
    0.5, 0.5, "synthetic incomplete Cartesian grid");
catch exception
  didRejectIncompleteGrid = ...
    exception.identifier == "audit_native_grid_pairs:IncompleteNativeGrid";
end
assert(didRejectIncompleteGrid, ...
  "The native-grid audit accepted a duplicate coordinate and missing site.");
generatorText = string(fileread(which( ...
  "generate_grain_boundary_misorientation_distribution")));
assert(contains(generatorText, "calcGrains(ebsdFull, 'unitCell'"));
forbiddenCalls = [
  "calcGrains(ebsd('indexed')";
  "calcGrains(ebsd(""indexed"")";
  "calcGrains(ebsdFull('indexed')";
  "calcGrains(ebsdFull(""indexed"")";
  "gridify(";
  "interpolate(";
  "smooth("
];
assert(~any(contains(generatorText, forbiddenCalls)), ...
  "The generator contains an obsolete reconstruction or filtering call.");
assert(contains(generatorText, "partition_ti_hex_boundary_segments("));
assert(contains(generatorText, "audit_native_grid_pairs("));

weights = [1;3;2];
edges = 1:0.5:94;
[stats, distribution] = summarize_weighted_boundary_angles( ...
  [1.2;2.2;15.2], weights, edges);
assert(stats.segment_count == 3);
assert(abs(stats.total_boundary_length_um - 6) < 1e-12);
assert(abs(stats.weighted_mean_angle_deg - 38.2/6) < 1e-12);
assert(abs(stats.weighted_median_angle_deg - 2.2) < 1e-12);
assert(abs(sum(distribution.length_fraction) - 1) < 1e-12);
assert(abs(sum(distribution.probability_density_per_degree .* ...
  distribution.bin_width_deg) - 1) < 1e-12);

if scanRoot ~= ""
  assert(outputDir ~= "", "outputDir is required for integration testing.");
  expectedFiles = [
    "grain_boundary_detection_sensitivity.csv";
    "grain_boundary_misorientation_distribution.csv";
    "grain_boundary_misorientation_distribution.png";
    "grain_boundary_misorientation_metrics.png";
    "grain_boundary_misorientation_summary.csv"
  ];
  prepare_empty_output_folder(outputDir, expectedFiles);
  generationStart = datetime("now");

  [summary, sensitivity, fullDistribution] = ...
    generate_grain_boundary_misorientation_distribution(scanRoot, outputDir);

  sensitivitySchema = [
    "sample"; "folder"; "input_file"; "cold_reduction_percent";
    "detection_floor_deg"; "classification_angle_deg"; "min_pixel";
    "boundary_network_estimator"; "eligible_boundary_sources";
    "scan_unit"; "native_x_cell_count"; "native_y_cell_count";
    "native_x_step_um"; "native_y_step_um"; "pixel_area_um2";
    "mapped_pixel_count"; "mapped_area_um2"; "unindexed_pixel_count";
    "indexed_ti_hex_pixel_count";
    "indexed_ti_hex_area_um2"; "indexed_ti_hex_fraction_of_mapped_area";
    "ti_hex_grain_count";
    "inner_ti_ti_segment_count"; "inner_ti_ti_boundary_length_um";
    "outer_ti_ti_segment_count"; "outer_ti_ti_boundary_length_um";
    "audit_inner_lagb_segment_count"; "audit_inner_lagb_length_um";
    "audit_inner_high_angle_segment_count";
    "audit_inner_high_angle_length_um";
    "audit_outer_low_angle_segment_count"; "audit_outer_low_angle_length_um";
    "audit_outer_hagb_segment_count"; "audit_outer_hagb_length_um";
    "source_network_segment_count"; "source_network_boundary_length_um";
    "excluded_cross_class_segment_count";
    "excluded_cross_class_boundary_length_um";
    "native_ti_ti_endpoint_pair_count";
    "nonlocal_ti_ti_endpoint_pair_count";
    "max_ti_ti_endpoint_distance_um"; "eligible_segment_count";
    "total_eligible_boundary_length_um"; "lagb_segment_count";
    "hagb_segment_count"; "lagb_length_um"; "hagb_length_um";
    "lagb_length_fraction"; "hagb_length_fraction";
    "total_eligible_length_density_um_per_um2";
    "lagb_length_density_um_per_um2"; "hagb_length_density_um_per_um2";
    "weighted_mean_angle_deg"; "weighted_median_angle_deg";
    "min_angle_deg"; "max_angle_deg"
  ];
  intervalSchema = [
    "length_1_2_um"; "fraction_1_2"; "length_2_5_um"; "fraction_2_5";
    "length_5_10_um"; "fraction_5_10"; "length_10_15_um";
    "fraction_10_15"; "length_15_94_um"; "fraction_15_94"
  ];
  distributionSchema = [
    "sample"; "folder"; "input_file"; "cold_reduction_percent";
    "detection_floor_deg"; "classification_angle_deg"; "min_pixel";
    "boundary_network_estimator"; "eligible_boundary_sources";
    "scan_unit"; "indexed_ti_hex_area_um2"; "bin_lower_deg";
    "bin_upper_deg"; "bin_width_deg"; "bin_center_deg"; "segment_count";
    "boundary_length_um"; "length_fraction";
    "probability_density_per_degree"; "cumulative_length_fraction"
  ];

  assert(height(summary) == 6);
  assert(height(sensitivity) == 18);
  assert(height(fullDistribution) == 1116);
  assert(isequal(string(sensitivity.Properties.VariableNames)', ...
    sensitivitySchema));
  assert(isequal(string(summary.Properties.VariableNames)', ...
    [sensitivitySchema; intervalSchema]));
  assert(isequal(string(fullDistribution.Properties.VariableNames)', ...
    distributionSchema));
  expectedSamples = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"];
  expectedInputFiles = [
    "ebsd_sample_7_map_15.ctf";
    "ebsd_sample_648_map_13.ctf";
    "ebsd_sample_602_map_11.ctf";
    "ebsd_sample_56_map_9.ctf";
    "ebsd_sample_525_map_7.ctf";
    "ebsd_sample_5_map_3.ctf"
  ];
  assert(isequal(summary.sample, expectedSamples));
  assert(isequal(summary.input_file, expectedInputFiles));
  assert(~any(contains(sensitivity.input_file, "_denoised")));
  assert(all(sensitivity.min_pixel == 5));
  assert(isequal(unique(sensitivity.detection_floor_deg), [0.5;1;2]));
  assert(all(sensitivity.boundary_network_estimator == ...
    "MTEX thresholded Ti-Hex/Ti-Hex boundary-segment network"));
  assert(all(sensitivity.eligible_boundary_sources == ...
    "innerBoundary LAGB + grains.boundary HAGB"));
  assert_scan_audit(sensitivity);
  assert_source_angle_conservation(sensitivity);

  primaryRows = sensitivity.detection_floor_deg == 1;
  assert(isequal(summary(:, cellstr(sensitivitySchema)), ...
    sensitivity(primaryRows, :)));
  probeRow = sensitivity.sample == "6.48d" & primaryRows;
  assert(nnz(probeRow) == 1);
  assert(sensitivity.ti_hex_grain_count(probeRow) == 2178);
  for sample = expectedSamples'
    sampleRows = sensitivity.sample == sample;
    assert(nnz(sampleRows) == 3);
    assert(isequal(sort(sensitivity.detection_floor_deg(sampleRows)), ...
      [0.5;1;2]));
  end

  intervalLength = summary.length_1_2_um + summary.length_2_5_um + ...
    summary.length_5_10_um + summary.length_10_15_um + ...
    summary.length_15_94_um;
  intervalFraction = summary.fraction_1_2 + summary.fraction_2_5 + ...
    summary.fraction_5_10 + summary.fraction_10_15 + ...
    summary.fraction_15_94;
  assert(all(abs(intervalLength - ...
    summary.total_eligible_boundary_length_um) < 1e-8));
  assert(all(abs(intervalFraction - 1) < 1e-10));
  assert(all(abs(summary.length_15_94_um - summary.hagb_length_um) < 1e-8));
  assert(all(abs(summary.length_1_2_um + summary.length_2_5_um + ...
    summary.length_5_10_um + summary.length_10_15_um - ...
    summary.lagb_length_um) < 1e-8));
  assert_distribution_conservation(summary, fullDistribution, expectedSamples);

  summaryRead = readtable(fullfile(outputDir, ...
    "grain_boundary_misorientation_summary.csv"), "TextType", "string");
  sensitivityRead = readtable(fullfile(outputDir, ...
    "grain_boundary_detection_sensitivity.csv"), "TextType", "string");
  distributionRead = readtable(fullfile(outputDir, ...
    "grain_boundary_misorientation_distribution.csv"), "TextType", "string");
  assert(isequal(string(summaryRead.Properties.VariableNames)', ...
    [sensitivitySchema; intervalSchema]));
  assert(isequal(string(sensitivityRead.Properties.VariableNames)', ...
    sensitivitySchema));
  assert(isequal(string(distributionRead.Properties.VariableNames)', ...
    distributionSchema));
  assert(isequal(summaryRead.sample, expectedSamples));
  assert_scan_audit(sensitivityRead);
  assert_source_angle_conservation(sensitivityRead);
  assert_distribution_conservation(summaryRead, distributionRead, ...
    expectedSamples);
  assert_fresh_exact_artifacts(outputDir, expectedFiles, generationStart);
end

fprintf("GRAIN_BOUNDARY_MISORIENTATION_TESTS_OK\n");
end

function prepare_empty_output_folder(outputDir, expectedFiles)
if ~isfolder(outputDir)
  mkdir(outputDir);
end
for fileName = expectedFiles'
  filePath = fullfile(outputDir, fileName);
  if isfile(filePath)
    delete(filePath);
  end
end
remaining = dir(outputDir);
remaining = remaining(~[remaining.isdir]);
assert(isempty(remaining), ...
  "Integration output folder must contain no stale or unrelated files.");
end

function assert_fresh_exact_artifacts(outputDir, expectedFiles, generationStart)
produced = dir(outputDir);
produced = produced(~[produced.isdir]);
producedNames = reshape(sort(string({produced.name})), [], 1);
assert(isequal(producedNames, sort(expectedFiles)), ...
  "The generator did not create the exact five-artifact output set.");
for fileIndex = 1:numel(expectedFiles)
  fileInfo = dir(fullfile(outputDir, expectedFiles(fileIndex)));
  assert(isscalar(fileInfo) && fileInfo.bytes > 100, ...
    "Generated artifact is missing or empty: %s", expectedFiles(fileIndex));
  fileTimestamp = datetime(fileInfo.datenum, "ConvertFrom", "datenum");
  assert(fileTimestamp >= generationStart - seconds(2), ...
    "Generated artifact is not fresh: %s", expectedFiles(fileIndex));
end
end

function assert_scan_audit(rows)
tol = 1e-10;
assert(all(rows.scan_unit == "um"));
assert(all(rows.native_x_cell_count == 600));
assert(all(rows.native_y_cell_count == 600));
assert(all(abs(rows.native_x_step_um - 0.5) < tol));
assert(all(abs(rows.native_y_step_um - 0.5) < tol));
assert(all(rows.mapped_pixel_count == 360000));
assert(all(rows.unindexed_pixel_count > 0 & ...
  rows.unindexed_pixel_count < rows.mapped_pixel_count));
assert(all(rows.indexed_ti_hex_pixel_count <= ...
  rows.mapped_pixel_count - rows.unindexed_pixel_count));
assert(all(abs(rows.pixel_area_um2 - 0.25) < tol));
assert(all(abs(rows.mapped_area_um2 - 90000) < tol));
assert(all(abs(rows.indexed_ti_hex_area_um2 - ...
  rows.indexed_ti_hex_pixel_count .* rows.pixel_area_um2) < 1e-8));
assert(all(abs(rows.indexed_ti_hex_fraction_of_mapped_area - ...
  rows.indexed_ti_hex_area_um2 ./ rows.mapped_area_um2) < tol));
assert(all(rows.ti_hex_grain_count > 0 & ...
  rows.ti_hex_grain_count == fix(rows.ti_hex_grain_count)));
assert(all(rows.ti_hex_grain_count <= rows.indexed_ti_hex_pixel_count));
assert(all(rows.min_angle_deg >= rows.detection_floor_deg - tol));
assert(all(rows.max_angle_deg <= 94 + tol));
assert(all(rows.native_ti_ti_endpoint_pair_count == ...
  rows.source_network_segment_count));
assert(all(rows.nonlocal_ti_ti_endpoint_pair_count == 0));
assert(all(rows.max_ti_ti_endpoint_distance_um <= ...
  max(rows.native_x_step_um, rows.native_y_step_um) + tol));
end

function assert_source_angle_conservation(rows)
tol = 1e-8;
assert(all(rows.inner_ti_ti_segment_count == ...
  rows.audit_inner_lagb_segment_count + ...
  rows.audit_inner_high_angle_segment_count));
assert(all(rows.outer_ti_ti_segment_count == ...
  rows.audit_outer_low_angle_segment_count + ...
  rows.audit_outer_hagb_segment_count));
assert(all(rows.source_network_segment_count == ...
  rows.inner_ti_ti_segment_count + rows.outer_ti_ti_segment_count));
assert(all(rows.eligible_segment_count == ...
  rows.audit_inner_lagb_segment_count + ...
  rows.audit_outer_hagb_segment_count));
assert(all(rows.excluded_cross_class_segment_count == ...
  rows.audit_inner_high_angle_segment_count + ...
  rows.audit_outer_low_angle_segment_count));
assert(all(rows.source_network_segment_count == rows.eligible_segment_count + ...
  rows.excluded_cross_class_segment_count));

assert(all(abs(rows.inner_ti_ti_boundary_length_um - ...
  rows.audit_inner_lagb_length_um - ...
  rows.audit_inner_high_angle_length_um) < tol));
assert(all(abs(rows.outer_ti_ti_boundary_length_um - ...
  rows.audit_outer_low_angle_length_um - ...
  rows.audit_outer_hagb_length_um) < tol));
assert(all(abs(rows.source_network_boundary_length_um - ...
  rows.inner_ti_ti_boundary_length_um - ...
  rows.outer_ti_ti_boundary_length_um) < tol));
assert(all(abs(rows.total_eligible_boundary_length_um - ...
  rows.audit_inner_lagb_length_um - rows.audit_outer_hagb_length_um) < tol));
assert(all(abs(rows.excluded_cross_class_boundary_length_um - ...
  rows.audit_inner_high_angle_length_um - ...
  rows.audit_outer_low_angle_length_um) < tol));
assert(all(abs(rows.source_network_boundary_length_um - ...
  rows.total_eligible_boundary_length_um - ...
  rows.excluded_cross_class_boundary_length_um) < tol));

assert(all(rows.lagb_segment_count == rows.audit_inner_lagb_segment_count));
assert(all(rows.hagb_segment_count == rows.audit_outer_hagb_segment_count));
assert(all(abs(rows.lagb_length_um - rows.audit_inner_lagb_length_um) < tol));
assert(all(abs(rows.hagb_length_um - rows.audit_outer_hagb_length_um) < tol));
assert(all(abs(rows.lagb_length_um + rows.hagb_length_um - ...
  rows.total_eligible_boundary_length_um) < tol));
assert(all(abs(rows.lagb_length_fraction + rows.hagb_length_fraction - 1) ...
  < 1e-10));
assert(all(abs(rows.total_eligible_length_density_um_per_um2 - ...
  rows.lagb_length_density_um_per_um2 - ...
  rows.hagb_length_density_um_per_um2) < 1e-12));
assert(all(abs(rows.total_eligible_length_density_um_per_um2 - ...
  rows.total_eligible_boundary_length_um ./ ...
  rows.indexed_ti_hex_area_um2) < 1e-12));
end

function assert_distribution_conservation(summary, distribution, sampleNames)
for sample = sampleNames'
  summaryRow = summary.sample == sample;
  distributionRows = distribution.sample == sample;
  assert(nnz(summaryRow) == 1);
  assert(nnz(distributionRows) == 186);
  assert(all(distribution.detection_floor_deg(distributionRows) == 1));
  assert(all(distribution.classification_angle_deg(distributionRows) == ...
    summary.classification_angle_deg(summaryRow)));
  assert(all(distribution.min_pixel(distributionRows) == ...
    summary.min_pixel(summaryRow)));
  assert(all(distribution.cold_reduction_percent(distributionRows) == ...
    summary.cold_reduction_percent(summaryRow)));
  assert(all(distribution.input_file(distributionRows) == ...
    summary.input_file(summaryRow)));
  assert(all(distribution.boundary_network_estimator(distributionRows) == ...
    summary.boundary_network_estimator(summaryRow)));
  assert(all(distribution.eligible_boundary_sources(distributionRows) == ...
    summary.eligible_boundary_sources(summaryRow)));
  assert(all(distribution.scan_unit(distributionRows) == ...
    summary.scan_unit(summaryRow)));
  assert(all(abs(distribution.indexed_ti_hex_area_um2(distributionRows) - ...
    summary.indexed_ti_hex_area_um2(summaryRow)) < 1e-10));
  assert(sum(distribution.segment_count(distributionRows)) == ...
    summary.eligible_segment_count(summaryRow));
  assert(abs(sum(distribution.boundary_length_um(distributionRows)) - ...
    summary.total_eligible_boundary_length_um(summaryRow)) < 1e-8);
  assert(abs(sum(distribution.length_fraction(distributionRows)) - 1) ...
    < 1e-10);
  assert(abs(sum(distribution.probability_density_per_degree( ...
    distributionRows) .* distribution.bin_width_deg(distributionRows)) ...
    - 1) < 1e-10);
  assert(abs(distribution.cumulative_length_fraction( ...
    find(distributionRows, 1, "last")) - 1) < 1e-10);
end
end
