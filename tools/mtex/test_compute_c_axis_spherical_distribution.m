function test_compute_c_axis_spherical_distribution()
%TEST_COMPUTE_C_AXIS_SPHERICAL_DISTRIBUTION Focused spherical c-axis tests.

test_equal_area_grid();
test_spherical_distribution();
fprintf('test_compute_c_axis_spherical_distribution passed.\n');
end

function test_equal_area_grid()
[grid,gridVectors] = build_c_axis_equal_area_grid(18,72);
assert(height(grid) == 18*72);
assert(isequal(string(grid.Properties.VariableNames), ...
  ["grid_index","theta_ad_deg","phi_about_ad_deg","cell_weight"]));
assert(length(gridVectors) == height(grid));
assert(abs(sum(grid.cell_weight)-1) < 1e-12);
assert(all(grid.theta_ad_deg > 0 & grid.theta_ad_deg < 90));
assert(all(grid.phi_about_ad_deg >= -180 & ...
  grid.phi_about_ad_deg < 180));
end

function test_spherical_distribution()
xyz = [1 0 0;1 0 0;-1 0 0];
weights = [1;2;3];
[parallel,~,audit] = compute_c_axis_spherical_distribution( ...
  xyz,weights,36,144,5);
assert(abs(sum(parallel.cell_weight .* parallel.mrd)-1) < 5e-3);
assert(audit.maximum_theta_ad_deg < 5);
assert(isfield(audit,"grid_maximum_mrd"));
assert(isfield(audit,"grid_maximum_theta_ad_deg"));
assert(isfield(audit,"grid_maximum_phi_about_ad_deg"));
assert(audit.grid_maximum_theta_ad_deg > 5);
assert(audit.canonicalized_source_count == 1);

[positive,~,~] = compute_c_axis_spherical_distribution( ...
  [0.4 0.8 0.2],[1],18,72,7.5);
[negative,~,~] = compute_c_axis_spherical_distribution( ...
  [-0.4 -0.8 -0.2],[1],18,72,7.5);
assert(max(abs(positive.mrd-negative.mrd)) < 1e-10);

[grid,gridVectors] = build_c_axis_equal_area_grid(36,144);
randomXYZ = [gridVectors.x(:),gridVectors.y(:),gridVectors.z(:)];
random = compute_c_axis_spherical_distribution( ...
  randomXYZ,grid.cell_weight,36,144,10);
assert(max(abs(random.mrd-1)) < 0.15);
end
