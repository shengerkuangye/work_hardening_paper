function test_grain_boundary_misorientation_distribution(scanRoot, outputDir)
arguments
  scanRoot (1,1) string = ""
  outputDir (1,1) string = ""
end

componentId = [1;1;2;3;3];
segLengthUm = [0.4;0.7;0.9;0.25;0.25];
[keepMask, componentLength] = component_length_mask( ...
  componentId, segLengthUm, 1.0);
assert(isequal(keepMask, logical([1;1;0;0;0])));
assert(isequal(componentLength, [1.1;0.9;0.5]));
[emptyKeep, emptyLength] = component_length_mask([], [], 1.0);
assert(isempty(emptyKeep) && isempty(emptyLength));

thetaDeg = [1.2;2.2;15.2];
weights = [1;3;2];
edges = 1:0.5:94;
[stats, distribution] = summarize_weighted_boundary_angles( ...
  thetaDeg, weights, edges);
assert(stats.segment_count == 3);
assert(abs(stats.total_boundary_length_um - 6) < 1e-12);
assert(abs(stats.weighted_mean_angle_deg - 38.2/6) < 1e-12);
assert(abs(stats.weighted_median_angle_deg - 2.2) < 1e-12);
assert(abs(sum(distribution.length_fraction) - 1) < 1e-12);
assert(abs(sum(distribution.probability_density_per_degree .* ...
  distribution.bin_width_deg) - 1) < 1e-12);

if scanRoot ~= ""
  assert(outputDir ~= "", "outputDir is required for integration testing.");
  [summary, sensitivity, fullDistribution] = ...
    generate_grain_boundary_misorientation_distribution(scanRoot, outputDir);
  assert(height(summary) == 6);
  assert(height(sensitivity) == 18);
  assert(height(fullDistribution) == 1116);
  expectedSamples = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"];
  assert(isequal(summary.sample, expectedSamples));
  assert(all(summary.total_retained_boundary_length_um > 0));
  assert(all(abs(summary.lagb_length_um + summary.hagb_length_um - ...
    summary.total_retained_boundary_length_um) < 1e-8));
  intervalSum = summary.length_1_2_um + summary.length_2_5_um + ...
    summary.length_5_10_um + summary.length_10_15_um + ...
    summary.length_15_94_um;
  assert(all(abs(intervalSum - summary.total_retained_boundary_length_um) ...
    < 1e-8));
  assert(isequal(unique(sensitivity.detection_floor_deg), [0.5;1;2]));
  for sample = expectedSamples'
    rows = fullDistribution.sample == sample;
    assert(nnz(rows) == 186);
    assert(abs(sum(fullDistribution.length_fraction(rows)) - 1) < 1e-10);
  end
  expectedFiles = [
    "grain_boundary_misorientation_distribution.csv";
    "grain_boundary_misorientation_summary.csv";
    "grain_boundary_detection_sensitivity.csv";
    "grain_boundary_misorientation_distribution.png";
    "grain_boundary_misorientation_metrics.png"
  ];
  assert(all(isfile(fullfile(outputDir, expectedFiles))));
end

fprintf("GRAIN_BOUNDARY_MISORIENTATION_TESTS_OK\n");
end
