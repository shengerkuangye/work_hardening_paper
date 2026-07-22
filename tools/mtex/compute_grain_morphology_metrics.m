function [byGrain, summary, quantiles] = ...
  compute_grain_morphology_metrics(grains, meta, options)
%COMPUTE_GRAIN_MORPHOLOGY_METRICS Unsmoothed MTEX grain-shape metrics.

arguments
  grains grain2d
  meta table
  options (1,1) struct
end
assert(height(meta) == 1 && all(ismember( ...
  ["sample","diameter_mm","cold_reduction_percent","variant"], ...
  string(meta.Properties.VariableNames))));
assert(isfield(options, "grain_detection_deg") && ...
  isfield(options, "min_grain_pixels"));
grainDetectionDeg = options.grain_detection_deg;
minGrainPixels = options.min_grain_pixels;
assert(isscalar(grainDetectionDeg) && isfinite(grainDetectionDeg) && ...
  grainDetectionDeg > 0);
assert(isscalar(minGrainPixels) && isfinite(minGrainPixels) && ...
  minGrainPixels >= 1 && minGrainPixels == fix(minGrainPixels));

keep = grains.numPixel >= minGrainPixels;
measured = grains(keep);
assert(~isempty(measured), "No grains satisfy the minimum-pixel rule.");

areaUm2 = double(area(measured));
areaUm2 = areaUm2(:);
perimeterUm = double(perimeter(measured));
perimeterUm = perimeterUm(:);
[~, ellipseLong, ellipseShort] = fitEllipse(measured);
ellipseLongUm = 2 * double(norm(ellipseLong(:)));
ellipseShortUm = 2 * double(norm(ellipseShort(:)));
aspectRatio = ellipseLongUm ./ ellipseShortUm;
longDirection = longAxis(measured);
longAngleDeg = fold_antipodal_angle_to_ad( ...
  atan2d(double(longDirection.y(:)), double(longDirection.x(:))));
maxFeretUm = double(norm(caliper(measured, 'longest')));
minFeretUm = double(norm(caliper(measured, 'shortest')));

grain_id = double(measured.id(:));
num_pixel = double(measured.numPixel(:));
boundary_touching = logical(isBoundary(measured));
boundary_touching = boundary_touching(:);
area_um2 = areaUm2;
ecd_um = 2 * sqrt(areaUm2 / pi);
perimeter_um = perimeterUm;
ellipse_long_axis_um = ellipseLongUm;
ellipse_short_axis_um = ellipseShortUm;
aspect_ratio = aspectRatio;
long_axis_ad_angle_deg = longAngleDeg;
max_feret_um = maxFeretUm;
min_feret_um = minFeretUm;
shape_factor = double(shapeFactor(measured));
shape_factor = shape_factor(:);
byGrain = table(grain_id, num_pixel, boundary_touching, area_um2, ...
  ecd_um, perimeter_um, ellipse_long_axis_um, ...
  ellipse_short_axis_um, aspect_ratio, long_axis_ad_angle_deg, ...
  max_feret_um, min_feret_um, shape_factor);
byGrain = add_common_columns(byGrain, meta);

summary = table(string(meta.sample), double(meta.diameter_mm), ...
  double(meta.cold_reduction_percent), string(meta.variant), ...
  grainDetectionDeg, minGrainPixels, height(byGrain), sum(areaUm2), ...
  median(ecd_um), weighted_quantile(ecd_um, areaUm2, 0.5), ...
  median(aspect_ratio), weighted_quantile(aspect_ratio, areaUm2, 0.5), ...
  median(max_feret_um), median(min_feret_um), median(shape_factor), ...
  weighted_quantile(long_axis_ad_angle_deg, areaUm2, 0.5), ...
  mean(boundary_touching), 'VariableNames', cellstr( ...
  comprehensive_ebsd_output_contract().summaryColumns. ...
  grain_morphology_summary));

quantiles = morphology_quantiles(byGrain);
end

function values = fold_antipodal_angle_to_ad(values)
values = mod(values, 180);
values = min(values, 180 - values);
angleToleranceDeg = 1e-10;
values(abs(values) < angleToleranceDeg) = 0;
values(abs(values - 90) < angleToleranceDeg) = 90;
assert(all(values >= 0 & values <= 90));
end

function output = add_common_columns(input, meta)
n = height(input);
output = addvars(input, repmat(string(meta.sample), n, 1), ...
  repmat(double(meta.diameter_mm), n, 1), ...
  repmat(double(meta.cold_reduction_percent), n, 1), ...
  repmat(string(meta.variant), n, 1), 'Before', 1, ...
  'NewVariableNames', {'sample','diameter_mm', ...
  'cold_reduction_percent','variant'});
end

function output = morphology_quantiles(byGrain)
probabilities = [0.1;0.25;0.5;0.75;0.9];
scopeNames = ["all";"exclude_boundary_touching"];
weightingNames = ["number";"area"];
metricNames = ["ecd_um","aspect_ratio","max_feret_um", ...
  "min_feret_um","shape_factor","long_axis_ad_angle_deg"];
rows = cell(0, 3 + numel(metricNames));
for scopeIndex = 1:numel(scopeNames)
  if scopeNames(scopeIndex) == "all"
    mask = true(height(byGrain), 1);
  else
    mask = ~byGrain.boundary_touching;
  end
  for weightingIndex = 1:numel(weightingNames)
    for probabilityIndex = 1:numel(probabilities)
      row = cell(1, 3 + numel(metricNames));
      row(1:3) = {scopeNames(scopeIndex), ...
        weightingNames(weightingIndex), probabilities(probabilityIndex)};
      for metricIndex = 1:numel(metricNames)
        values = byGrain.(metricNames(metricIndex))(mask);
        if weightingNames(weightingIndex) == "number"
          row{3 + metricIndex} = unweighted_quantile( ...
            values, probabilities(probabilityIndex));
        else
          row{3 + metricIndex} = weighted_quantile(values, ...
            byGrain.area_um2(mask), probabilities(probabilityIndex));
        end
      end
      rows(end + 1, :) = row; %#ok<AGROW>
    end
  end
end
output = cell2table(rows, 'VariableNames', ...
  cellstr(["scope","weighting","quantile_probability",metricNames]));
output.scope = string(output.scope);
output.weighting = string(output.weighting);
numericNames = ["quantile_probability", metricNames];
for name = numericNames
  if iscell(output.(name))
    output.(name) = cell2mat(output.(name));
  end
end
end

function q = unweighted_quantile(values, probability)
values = double(values(:));
values = values(isfinite(values));
if isempty(values)
  q = NaN;
else
  q = quantile(values, probability);
end
end

function q = weighted_quantile(values, weights, probability)
values = double(values(:));
weights = double(weights(:));
valid = isfinite(values) & isfinite(weights) & weights > 0;
values = values(valid);
weights = weights(valid);
if isempty(values)
  q = NaN;
  return
end
[values, order] = sort(values);
weights = weights(order);
index = find(cumsum(weights) >= probability * sum(weights), 1);
q = values(index);
end
