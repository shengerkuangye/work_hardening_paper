function test_neighbor_misorientation_distribution(scanRoot, outputDir)
arguments
  scanRoot (1,1) string = ""
  outputDir (1,1) string = ""
end

phaseMask = logical([1 1 0; 1 1 1]);
[index1, index2] = four_neighbor_pairs(phaseMask);
pairs = sort([index1, index2], 2);
expectedPairs = sort([1 3; 2 4; 4 6; 1 2; 3 4], 2);
assert(size(pairs,1) == 5);
assert(size(unique(pairs, "rows"),1) == 5);
assert(isequal(sortrows(pairs), sortrows(expectedPairs)));

binEdgesDeg = 0:0.1:94;
thetaDeg = [0.04; 0.14; 1.04; 93.84];
[stats, distribution] = summarize_misorientation_angles( ...
  thetaDeg, binEdgesDeg);
assert(stats.pair_count == 4);
assert(abs(stats.mean_angle_deg - mean(thetaDeg)) < 1e-12);
assert(abs(stats.mode_bin_center_deg - 0.05) < 1e-12);
assert(sum(distribution.pair_count) == 4);
assert(abs(sum(distribution.probability) - 1) < 1e-12);
assert(abs(distribution.cumulative_probability(end) - 1) < 1e-12);
assert(abs(sum(distribution.probability_density_per_degree .* ...
  distribution.bin_width_deg) - 1) < 1e-12);

if scanRoot ~= ""
  assert(outputDir ~= "", "outputDir is required for integration testing.");
  [summary, fullDistribution] = ...
    generate_neighbor_misorientation_distribution(scanRoot, outputDir);
  assert(height(summary) == 7);
  assert(height(fullDistribution) == 6580);
  assert(isequal(summary.sample, ...
    ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d";"pooled"]));
  assert(all(summary.pair_count(1:6) > 0 & ...
    summary.pair_count(1:6) <= 718800));
  assert(summary.pair_count(7) == sum(summary.pair_count(1:6)));
  for sample = summary.sample'
    rows = fullDistribution.sample == sample;
    assert(nnz(rows) == 940);
    assert(abs(sum(fullDistribution.probability(rows)) - 1) < 1e-10);
    assert(abs(fullDistribution.cumulative_probability(find(rows,1,"last")) ...
      - 1) < 1e-10);
  end
  assert(isfile(fullfile(outputDir, ...
    "neighbor_misorientation_distribution.csv")));
  assert(isfile(fullfile(outputDir, ...
    "neighbor_misorientation_summary.csv")));
  assert(isfile(fullfile(outputDir, ...
    "neighbor_misorientation_distribution.png")));
end

fprintf("NEIGHBOR_MISORIENTATION_TESTS_OK\n");
end
