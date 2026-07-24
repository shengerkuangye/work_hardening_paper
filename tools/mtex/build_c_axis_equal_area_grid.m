function [gridTable,gridVectors] = build_c_axis_equal_area_grid(nMu,nPhi)
%BUILD_C_AXIS_EQUAL_AREA_GRID Equal-area AD-positive hemisphere grid.

assert(isscalar(nMu) && nMu >= 4 && nMu == fix(nMu));
assert(isscalar(nPhi) && nPhi >= 8 && nPhi == fix(nPhi));
mu = ((1:nMu)-0.5)/nMu;
phiDeg = -180 + ((1:nPhi)-0.5) * (360/nPhi);
[phiGridDeg,muGrid] = meshgrid(phiDeg,mu);
thetaGridDeg = acosd(muGrid);
radial = sqrt(max(0,1-muGrid.^2));
x = muGrid;
y = radial .* cosd(phiGridDeg);
z = radial .* sind(phiGridDeg);
gridVectors = vector3d(x(:),y(:),z(:));
grid_index = (1:numel(x))';
theta_ad_deg = thetaGridDeg(:);
phi_about_ad_deg = phiGridDeg(:);
cell_weight = repmat(1/numel(x),numel(x),1);
gridTable = table(grid_index,theta_ad_deg,phi_about_ad_deg, ...
  cell_weight);
end
