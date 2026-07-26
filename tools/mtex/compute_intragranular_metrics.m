function [summary, byPixel, byGrain] = compute_intragranular_metrics( ...
  ebsdFull, grains, sampleMeta, options)
%COMPUTE_INTRAGRANULAR_METRICS Calculate KAM, GROD, and GOS proxies.
% KAM, GROD, and GOS are orientation-gradient/spread proxies; none is a
% direct measurement of dislocation density.

arguments
  ebsdFull EBSD
  grains grain2d
  sampleMeta (1,1) struct
  options (1,1) struct
end

requiredMeta = ["sample", "diameter_mm", "cold_reduction_percent", ...
  "variant"];
requiredOptions = ["grain_detection_deg", "min_grain_pixels", ...
  "kam_orders", "kam_thresholds_deg", "axis_min_grod_deg"];
assert(all(isfield(sampleMeta, requiredMeta)), ...
  "sampleMeta is missing a required field.");
assert(all(isfield(options, requiredOptions)), ...
  "options is missing a required field.");
assert(all(ismember(options.kam_orders, [1 2])) && ...
  all(options.kam_thresholds_deg > 0));
assert(options.axis_min_grod_deg > 0);
assert(isfield(ebsdFull.prop, "grainId"), ...
  "ebsdFull must contain grainId from unsmoothed reconstruction.");

if isa(ebsdFull, "EBSDsquare") || isa(ebsdFull, "EBSDhex")
  ebsdGrid = ebsdFull;
else
  ebsdGrid = ebsdFull.gridify;
end
assert(isfield(ebsdGrid.prop, "grainId"), ...
  "gridify did not preserve grainId.");

tiGrid = ebsdGrid("Ti-Hex");
assert(~isempty(tiGrid), "Ti-Hex phase is absent.");
tiPhaseId = unique(double(tiGrid.phaseId));
assert(isscalar(tiPhaseId));
tiMaskGrid = reshape(double(ebsdGrid.phaseId), size(ebsdGrid)) == tiPhaseId;

grod = calcGROD(ebsdFull, grains);
grodAngleDegAll = reshape(double(grod.angle / degree), [], 1);
tiFull = ebsdFull("Ti-Hex");
tiPhaseIdFull = unique(double(tiFull.phaseId));
tiMaskFull = reshape(double(ebsdFull.phaseId), [], 1) == tiPhaseIdFull;
grodValid = tiMaskFull & isfinite(grodAngleDegAll);
assert(any(grodValid), "No valid Ti-Hex GROD pixels.");

axisCrystal = grod.axis("noSymmetry");
axisSpecimen = ebsdFull.rotations .* axisCrystal;
[crystalX, crystalY, crystalZ] = vector_components(axisCrystal);
[specimenX, specimenY, specimenZ] = vector_components(axisSpecimen);
axisValid = grodValid & grodAngleDegAll >= options.axis_min_grod_deg & ...
  all(isfinite([crystalX crystalY crystalZ specimenX specimenY specimenZ]), 2);

pixelRows = find(tiMaskFull);
byPixel = table(double(ebsdFull.id(pixelRows)), ...
  double(ebsdFull.x(pixelRows)), double(ebsdFull.y(pixelRows)), ...
  double(ebsdFull.grainId(pixelRows)), grodAngleDegAll(pixelRows), ...
  axisValid(pixelRows), crystalX(pixelRows), crystalY(pixelRows), ...
  crystalZ(pixelRows), specimenX(pixelRows), specimenY(pixelRows), ...
  specimenZ(pixelRows), 'VariableNames', ...
  ["ebsd_id", "x_um", "y_um", "grain_id", "grod_angle_deg", ...
  "grod_axis_valid", "grod_axis_crystal_x", "grod_axis_crystal_y", ...
  "grod_axis_crystal_z", "grod_axis_specimen_x_ad", ...
  "grod_axis_specimen_y_td_rd", "grod_axis_specimen_z_nd"]);
byPixel{~byPixel.grod_axis_valid, 7:12} = NaN;

nRows = numel(options.kam_orders) * numel(options.kam_thresholds_deg);
summary = empty_summary(nRows);
row = 0;
for order = options.kam_orders(:)'
  for thresholdDeg = options.kam_thresholds_deg(:)'
    row = row + 1;
    kam = ebsdGrid.KAM("order", order, ...
      "threshold", thresholdDeg * degree);
    kamDegGrid = double(kam / degree);
    kamValidGrid = tiMaskGrid & isfinite(kamDegGrid);
    kamValues = kamDegGrid(kamValidGrid);
    assert(~isempty(kamValues), ...
      "No valid KAM pixels for order %d and threshold %g degree.", ...
      order, thresholdDeg);

    if isfield(ebsdGrid.prop, "oldId")
      gridSourceIds = double(ebsdGrid.oldId(:));
    else
      gridSourceIds = double(ebsdGrid.id(:));
    end
    [found, pixelPosition] = ismember( ...
      double(byPixel.ebsd_id), gridSourceIds);
    assert(all(found), ...
      "Native-grid source IDs do not match source EBSD IDs.");
    kamByPixel = kamDegGrid(pixelPosition);
    variableName = sprintf("kam_order%d_threshold%gdeg", order, thresholdDeg);
    variableName = matlab.lang.makeValidName(variableName);
    byPixel.(variableName) = kamByPixel;

    summary(row,:) = summary_row(sampleMeta, options, order, thresholdDeg, ...
      kamValues, grodAngleDegAll(grodValid), grains);
  end
end

byGrain = grain_summary(grains, ebsdFull, grodAngleDegAll, axisValid, ...
  crystalX, crystalY, crystalZ, specimenX, specimenY, specimenZ, options);
end

function summary = empty_summary(n)
sample = strings(n,1);
diameter_mm = zeros(n,1);
cold_reduction_percent = zeros(n,1);
variant = strings(n,1);
grain_detection_deg = zeros(n,1);
min_grain_pixels = zeros(n,1);
kam_order = zeros(n,1);
kam_threshold_deg = zeros(n,1);
valid_pixel_count = zeros(n,1);
kam_mean_deg = zeros(n,1);
kam_median_deg = zeros(n,1);
kam_p90_deg = zeros(n,1);
grod_mean_deg = zeros(n,1);
grod_median_deg = zeros(n,1);
grod_p90_deg = zeros(n,1);
gos_number_mean_deg = zeros(n,1);
gos_area_weighted_mean_deg = zeros(n,1);
summary = table(sample, diameter_mm, cold_reduction_percent, variant, ...
  grain_detection_deg, min_grain_pixels, kam_order, kam_threshold_deg, ...
  valid_pixel_count, kam_mean_deg, kam_median_deg, kam_p90_deg, ...
  grod_mean_deg, grod_median_deg, grod_p90_deg, ...
  gos_number_mean_deg, gos_area_weighted_mean_deg);
end

function row = summary_row(meta, options, order, thresholdDeg, ...
  kamValues, grodValues, grains)
tiGrains = grains("Ti-Hex");
useGrain = tiGrains.numPixel >= options.min_grain_pixels;
tiGrains = tiGrains(useGrain);
assert(~isempty(tiGrains), "No Ti-Hex grains pass min_grain_pixels.");
gosDeg = double(tiGrains.GOS(:) / degree);
areas = abs(double(area(tiGrains)));
areas = areas(:);
assert(all(isfinite(gosDeg)) && all(isfinite(areas) & areas > 0));

row = table(string(meta.sample), double(meta.diameter_mm), ...
  double(meta.cold_reduction_percent), string(meta.variant), ...
  double(options.grain_detection_deg), double(options.min_grain_pixels), ...
  double(order), double(thresholdDeg), numel(kamValues), mean(kamValues), ...
  median(kamValues), prctile(kamValues, 90), mean(grodValues), ...
  median(grodValues), prctile(grodValues, 90), mean(gosDeg), ...
  sum(gosDeg .* areas) / sum(areas), 'VariableNames', ...
  ["sample", "diameter_mm", "cold_reduction_percent", "variant", ...
  "grain_detection_deg", "min_grain_pixels", "kam_order", ...
  "kam_threshold_deg", "valid_pixel_count", "kam_mean_deg", ...
  "kam_median_deg", "kam_p90_deg", "grod_mean_deg", ...
  "grod_median_deg", "grod_p90_deg", "gos_number_mean_deg", ...
  "gos_area_weighted_mean_deg"]);
end

function result = grain_summary(grains, ebsd, grodDeg, axisValid, ...
  cx, cy, cz, sx, sy, sz, options)
tiGrains = grains("Ti-Hex");
tiGrains = tiGrains(tiGrains.numPixel >= options.min_grain_pixels);
n = length(tiGrains);
grain_id = double(tiGrains.id(:));
num_pixels = double(tiGrains.numPixel(:));
area_um2 = abs(double(area(tiGrains)));
area_um2 = area_um2(:);
gos_deg = double(tiGrains.GOS(:) / degree);
valid_grod_pixel_count = zeros(n,1);
grod_mean_deg = zeros(n,1);
grod_median_deg = zeros(n,1);
grod_p90_deg = zeros(n,1);
valid_grod_axis_pixel_count = zeros(n,1);
axisValues = NaN(n,8);

pixelGrainIds = double(ebsd.grainId(:));
for k = 1:n
  inGrain = pixelGrainIds == grain_id(k) & isfinite(grodDeg);
  values = grodDeg(inGrain);
  valid_grod_pixel_count(k) = numel(values);
  grod_mean_deg(k) = mean(values);
  grod_median_deg(k) = median(values);
  grod_p90_deg(k) = prctile(values,90);
  stable = inGrain & axisValid;
  valid_grod_axis_pixel_count(k) = nnz(stable);
  if nnz(stable) >= 3
    [axisValues(k,1:3), axisValues(k,4)] = ...
      axial_principal([cx(stable), cy(stable), cz(stable)]);
    [axisValues(k,5:7), axisValues(k,8)] = ...
      axial_principal([sx(stable), sy(stable), sz(stable)]);
  end
end

result = table(grain_id, num_pixels, area_um2, gos_deg, ...
  valid_grod_pixel_count, grod_mean_deg, grod_median_deg, grod_p90_deg, ...
  valid_grod_axis_pixel_count, axisValues(:,1), axisValues(:,2), ...
  axisValues(:,3), axisValues(:,4), axisValues(:,5), axisValues(:,6), ...
  axisValues(:,7), axisValues(:,8), 'VariableNames', ...
  ["grain_id", "num_pixels", "area_um2", "gos_deg", ...
  "valid_grod_pixel_count", "grod_mean_deg", "grod_median_deg", ...
  "grod_p90_deg", "valid_grod_axis_pixel_count", ...
  "grod_crystal_axis_principal_x", "grod_crystal_axis_principal_y", ...
  "grod_crystal_axis_principal_z", "grod_crystal_axis_concentration", ...
  "grod_specimen_axis_principal_x_ad", ...
  "grod_specimen_axis_principal_y_td_rd", ...
  "grod_specimen_axis_principal_z_nd", ...
  "grod_specimen_axis_concentration"]);
end

function [direction, concentration] = axial_principal(xyz)
tensor = (xyz' * xyz) / size(xyz,1);
[vectors, values] = eig(tensor, "vector");
[lambda, index] = max(real(values));
direction = real(vectors(:,index))';
if direction(1) < 0 || (direction(1) == 0 && direction(2) < 0)
  direction = -direction;
end
concentration = (3 * lambda - 1) / 2;
end

function [x, y, z] = vector_components(v)
x = double(v.x(:));
y = double(v.y(:));
z = double(v.z(:));
end
