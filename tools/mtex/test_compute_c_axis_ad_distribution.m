function test_compute_c_axis_ad_distribution()
%TEST_COMPUTE_C_AXIS_AD_DISTRIBUTION Verify axial PDF and MRD outputs.

edges = 0:2:90;
[parallel,audit] = compute_c_axis_ad_distribution( ...
  [0;0],[1;3],edges);
assert(abs(sum(parallel.observed_probability)-1) < 1e-12);
assert(parallel.observed_probability(1) == 1);
assert(all(parallel.observed_probability(2:end) == 0));
assert(abs(audit.valid_source_weight-4) < 1e-12);
assert(abs(audit.random_probability_sum-1) < 1e-12);
expectedRandomProbability = ...
  cosd(edges(1:end-1)')-cosd(edges(2:end)');
assert(max(abs(parallel.random_probability- ...
  expectedRandomProbability)) < 1e-15);

normal = compute_c_axis_ad_distribution([90;90],[1;1],edges);
assert(normal.observed_probability(end) == 1);

n = 200000;
mu = ((1:n)'-0.5)/n;
theta = acosd(mu);
random = compute_c_axis_ad_distribution(theta,ones(n,1),edges);
assert(max(abs(random.mrd-1)) < 0.02);

assert_error(@() compute_c_axis_ad_distribution( ...
  [-0.1;20],[1;1],edges));
assert_error(@() compute_c_axis_ad_distribution( ...
  [10;20],[1;0],edges));
assert_error(@() compute_c_axis_ad_distribution( ...
  [10;20],[1;1],0:2:88));

fprintf("COMPUTE_C_AXIS_AD_DISTRIBUTION_TESTS_OK\n");
end

function assert_error(functionHandle)
didError = false;
try
  functionHandle();
catch
  didError = true;
end
assert(didError);
end
