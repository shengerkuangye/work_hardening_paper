function test_c_axis_ad_alignment_metrics(scanRoot, outputDir)
%TEST_C_AXIS_AD_ALIGNMENT_METRICS Verify c-axis alignment calculations.

arguments
  scanRoot (1,1) string = ""
  outputDir (1,1) string = ""
end

edges = 0:2:90;

scaledAD = vector3d(0.2, 0, 0);
scaledNormal = vector3d(0, 0, 0.3);
assert(abs(c_axis_angles_to_ad(scaledAD)) < 1e-12);
assert(abs(c_axis_angles_to_ad(scaledNormal) - 90) < 1e-12);

[parallelStats, parallelDist] = compute_c_axis_ad_statistics( ...
  zeros(100,1), 100, edges);
assert(abs(parallelStats.alignment_factor - 1) < 1e-12);
assert(parallelStats.fraction_within_15deg == 1);
assert(sum(parallelDist.pixel_count) == 100);
assert(abs(sum(parallelDist.probability) - 1) < 1e-12);

[normalStats, normalDist] = compute_c_axis_ad_statistics( ...
  90 * ones(100,1), 100, edges);
assert(abs(normalStats.alignment_factor + 0.5) < 1e-12);
assert(normalStats.fraction_within_45deg == 0);
assert(sum(normalDist.pixel_count) == 100);

[mixedStats, ~] = compute_c_axis_ad_statistics([0; 90], 4, edges);
assert(abs(mixedStats.mean_cos2 - 0.5) < 1e-12);
assert(abs(mixedStats.alignment_factor - 0.25) < 1e-12);
assert(abs(mixedStats.alpha_ti_pixel_fraction - 0.5) < 1e-12);

if scanRoot ~= ""
  assert(outputDir ~= "", "outputDir is required for integration testing.");
  [summary, distribution] = generate_c_axis_ad_alignment_metrics( ...
    scanRoot, outputDir);
  assert(height(summary) == 6);
  assert(height(distribution) == 270);
  assert(all(summary.total_pixel_count == 360000));
  assert(all(summary.alpha_ti_pixel_count > 0));
  assert(all(summary.alignment_factor >= -0.5 & ...
    summary.alignment_factor <= 1));
  assert(all(summary.fraction_within_15deg >= 0 & ...
    summary.fraction_within_15deg <= 1));
  for sample = unique(distribution.sample, "stable")'
    rows = distribution.sample == sample;
    assert(abs(sum(distribution.probability(rows)) - 1) < 1e-10);
    assert(abs(sum(distribution.probability_density_per_degree(rows) .* ...
      distribution.bin_width_deg(rows)) - 1) < 1e-10);
  end
  assert(isfile(fullfile(outputDir, "c_axis_ad_alignment_summary.csv")));
  assert(isfile(fullfile(outputDir, "c_axis_ad_angle_distribution.csv")));
end

fprintf("C_AXIS_AD_ALIGNMENT_TESTS_OK\n");
end
