function test_c_axis_distribution_functions(textureDir)
%TEST_C_AXIS_DISTRIBUTION_FUNCTIONS Verify c-axis distribution outputs.

test_compute_c_axis_ad_distribution();
if nargin == 1
  assert(isfolder(string(textureDir)));
end
fprintf("C_AXIS_DISTRIBUTION_FUNCTIONS_TESTS_OK\n");
end
