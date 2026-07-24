function [distribution,densityFunction,audit] = ...
  compute_c_axis_spherical_distribution( ...
  cAxisXYZ,weights,nMu,nPhi,kernelHalfwidthDeg)
%COMPUTE_C_AXIS_SPHERICAL_DISTRIBUTION Antipodal c-axis MRD.

cAxisXYZ = double(cAxisXYZ);
weights = double(weights(:));
assert(size(cAxisXYZ,2) == 3 && size(cAxisXYZ,1) == numel(weights));
assert(all(isfinite(cAxisXYZ),"all"));
assert(all(isfinite(weights) & weights > 0));
assert(isscalar(kernelHalfwidthDeg) && kernelHalfwidthDeg > 0);
norms = vecnorm(cAxisXYZ,2,2);
assert(all(norms > 0));
cAxisXYZ = cAxisXYZ ./ norms;
tol = 1e-12;
flipMask = cAxisXYZ(:,1) < -tol | ...
  (abs(cAxisXYZ(:,1)) <= tol & cAxisXYZ(:,2) < -tol) | ...
  (abs(cAxisXYZ(:,1)) <= tol & abs(cAxisXYZ(:,2)) <= tol & ...
   cAxisXYZ(:,3) < 0);
cAxisXYZ(flipMask,:) = -cAxisXYZ(flipMask,:);

sourceVectors = vector3d(cAxisXYZ(:,1),cAxisXYZ(:,2),cAxisXYZ(:,3));
[distribution,gridVectors] = build_c_axis_equal_area_grid(nMu,nPhi);
kernel = S2DeLaValleePoussinKernel( ...
  "halfwidth",kernelHalfwidthDeg*degree);
densityFunction = calcDensity(sourceVectors,"weights",weights, ...
  "kernel",kernel,"antipodal");
rawGridMrd = double(densityFunction.eval(gridVectors));
rawGridMean = sum(distribution.cell_weight .* rawGridMrd);
assert(isfinite(rawGridMean) && rawGridMean > 0);
rawDensityFunction = densityFunction / rawGridMean;
rawGridMrd = rawGridMrd / rawGridMean;
rawMinimumMrd = min(rawGridMrd);
rawNegative = max(-rawGridMrd,0);
rawNegativeMass = sum(distribution.cell_weight .* rawNegative);
rawNegativeCellFraction = mean(rawGridMrd < 0);
maximumNegativeMass = 0.01;
assert(rawNegativeMass <= maximumNegativeMass, ...
  "Negative spherical-density mass %.17g exceeds limit %.17g.", ...
  rawNegativeMass,maximumNegativeMass);
gridMrd = max(rawGridMrd,0);
projectedGridMean = sum(distribution.cell_weight .* gridMrd);
assert(isfinite(projectedGridMean) && projectedGridMean > 0);
nonnegativeProjectionScale = 1/projectedGridMean;
gridMrd = gridMrd * nonnegativeProjectionScale;
densityFunction = S2FunHandle( ...
  @(v) max(double(rawDensityFunction.eval(v)),0) * ...
  nonnegativeProjectionScale,"antipodal");
distribution.mrd = gridMrd;
[gridMaximumMrd,gridMaximumIndex] = max(gridMrd);
[maximumMrd,maximumVector] = max(densityFunction);
maximumXYZ = [maximumVector.x,maximumVector.y,maximumVector.z];
maximumXYZ = maximumXYZ / norm(maximumXYZ);
maximumFlip = maximumXYZ(1) < -tol | ...
  (abs(maximumXYZ(1)) <= tol && maximumXYZ(2) < -tol) | ...
  (abs(maximumXYZ(1)) <= tol && abs(maximumXYZ(2)) <= tol && ...
   maximumXYZ(3) < 0);
if maximumFlip
  maximumXYZ = -maximumXYZ;
end

audit = struct();
audit.valid_source_count = size(cAxisXYZ,1);
audit.valid_source_weight = sum(weights);
audit.canonicalized_source_count = nnz(flipMask);
audit.grid_mean_mrd = sum(distribution.cell_weight .* gridMrd);
audit.raw_minimum_mrd = rawMinimumMrd;
audit.raw_negative_mass = rawNegativeMass;
audit.raw_negative_cell_fraction = rawNegativeCellFraction;
audit.maximum_negative_mass = maximumNegativeMass;
audit.nonnegative_projection_scale = nonnegativeProjectionScale;
audit.maximum_mrd = maximumMrd;
audit.maximum_theta_ad_deg = acosd(maximumXYZ(1));
audit.maximum_phi_about_ad_deg = atan2d(maximumXYZ(3),maximumXYZ(2));
audit.grid_maximum_mrd = gridMaximumMrd;
audit.grid_maximum_theta_ad_deg = ...
  distribution.theta_ad_deg(gridMaximumIndex);
audit.grid_maximum_phi_about_ad_deg = ...
  distribution.phi_about_ad_deg(gridMaximumIndex);
end
