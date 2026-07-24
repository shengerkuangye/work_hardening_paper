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

theta = [12;20;28];
phi = [0;55;120];
sharpXyz = [cosd(theta), ...
  sind(theta).*cosd(phi),sind(theta).*sind(phi)];
[sharp,sharpFunction,sharpAudit] = ...
  compute_c_axis_spherical_distribution( ...
  sharpXyz,[1;2;1],36,144,5);
assert(all(sharp.mrd >= 0));
assert(abs(sum(sharp.cell_weight .* sharp.mrd)-1) < 5e-12);
[~,sharpGridVectors] = build_c_axis_equal_area_grid(36,144);
sharpFunctionMrd = double(sharpFunction.eval(sharpGridVectors));
assert(all(sharpFunctionMrd >= 0));
assert(max(abs(sharpFunctionMrd-sharp.mrd)) < 0.01);
assert(abs(sum(sharp.cell_weight .* sharpFunctionMrd)-1) < 1e-5);
assert(sharpAudit.raw_minimum_mrd < 0);
assert(sharpAudit.raw_negative_mass > 0);
assert(sharpAudit.raw_negative_mass < 0.005);
assert(sharpAudit.raw_negative_cell_fraction > 0);
assert(sharpAudit.maximum_negative_mass == 0.01);
assert(sharpAudit.nonnegative_projection_scale < 1);

[grid,gridVectors] = build_c_axis_equal_area_grid(36,144);
randomXYZ = [gridVectors.x(:),gridVectors.y(:),gridVectors.z(:)];
random = compute_c_axis_spherical_distribution( ...
  randomXYZ,grid.cell_weight,36,144,10);
assert(max(abs(random.mrd-1)) < 0.15);
end
