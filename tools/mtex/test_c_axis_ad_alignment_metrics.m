function test_c_axis_ad_alignment_metrics
%TEST_C_AXIS_AD_ALIGNMENT_METRICS Verify c-axis alignment calculations.

edges = 0:2:90;

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

fprintf("C_AXIS_AD_ALIGNMENT_TESTS_OK\n");
end
