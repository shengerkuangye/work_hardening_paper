function test_grain_boundary_threshold_comparison(scanRoot, outputDir)
arguments
  scanRoot (1,1) string = ""
  outputDir (1,1) string = ""
end

thetaDeg = [2.1;4.9;5;10;14.9;15;30];
segLengthUm = [1;2;3;4;5;6;7];
floorDeg = 2;
indexedAreaUm2 = 100;

stats5 = calculate_boundary_threshold_metrics(thetaDeg, segLengthUm, ...
  floorDeg, 5, indexedAreaUm2);
stats15 = calculate_boundary_threshold_metrics(thetaDeg, segLengthUm, ...
  floorDeg, 15, indexedAreaUm2);

assert(stats5.total_eligible_boundary_length_um == 28);
assert(stats15.total_eligible_boundary_length_um == 28);
assert(stats5.lagb_length_um == 3);
assert(stats15.lagb_length_um == 15);
assert(stats5.hagb_length_um == 25);
assert(stats15.hagb_length_um == 13);
assert(abs(stats5.lagb_length_fraction + stats5.hagb_length_fraction - 1) ...
  < 1e-12);
assert(abs(stats15.lagb_length_fraction + stats15.hagb_length_fraction - 1) ...
  < 1e-12);
assert(stats5.lagb_length_um <= stats15.lagb_length_um);
assert(abs(stats5.lagb_length_density_um_per_um2 - 0.03) < 1e-12);
assert(abs(stats15.lagb_length_density_um_per_um2 - 0.15) < 1e-12);

if scanRoot ~= "" || outputDir ~= ""
  error("Integration comparison is not implemented yet.");
end

fprintf("GRAIN_BOUNDARY_THRESHOLD_COMPARISON_TESTS_OK\n");
end
