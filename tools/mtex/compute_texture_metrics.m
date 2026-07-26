function [summary, distribution, model] = compute_texture_metrics( ...
  pixelOrientations, grainMeanOrientations, grainAreas, sampleMeta, options)
%COMPUTE_TEXTURE_METRICS Calculate separately weighted alpha-Ti textures.

arguments
  pixelOrientations orientation
  grainMeanOrientations orientation
  grainAreas (:,1) double
  sampleMeta (1,1) struct
  options (1,1) struct
end

requiredMeta = ["sample", "diameter_mm", "cold_reduction_percent", ...
  "variant"];
requiredOptions = ["kernel_halfwidth_deg", "grid_resolution_deg", ...
  "c_axis_component_zero_tolerance", "min_grain_pixels"];
assert(all(isfield(sampleMeta, requiredMeta)));
assert(all(isfield(options, requiredOptions)));
assert(~isempty(pixelOrientations) && ~isempty(grainMeanOrientations));
grainAreas = grainAreas(:);
assert(numel(grainAreas) == length(grainMeanOrientations));
assert(all(isfinite(grainAreas) & grainAreas > 0));
assert(pixelOrientations.CS == grainMeanOrientations.CS);
assert(options.kernel_halfwidth_deg > 0 && options.grid_resolution_deg > 0);
assert(options.c_axis_component_zero_tolerance > 0);
assert(options.min_grain_pixels >= 1 && ...
  options.min_grain_pixels == fix(options.min_grain_pixels));

pixelIds = source_ids(options, "pixel_ids", length(pixelOrientations));
grainIds = source_ids(options, "grain_ids", length(grainMeanOrientations));
kernel = SO3DeLaValleePoussinKernel("halfwidth", ...
  options.kernel_halfwidth_deg * degree);
pixelRbfOdf = calcDensity(pixelOrientations, "kernel", kernel, ...
  "weights", ones(length(pixelOrientations),1), "silent");
% MTEX's weighted density uses the supplied weights on an absolute scale.
% Rescale to unit mean so the ODF integral remains one while preserving the
% relative grain-area weighting.
densityGrainWeights = grainAreas / mean(grainAreas);
grainRbfOdf = calcDensity(grainMeanOrientations, "kernel", kernel, ...
  "weights", densityGrainWeights, "silent");
harmonicBandwidth = kernel.bandwidth;
[pixelOdf, pixelDensityNormalization] = ...
  as_harmonic_odf(pixelRbfOdf, harmonicBandwidth);
[grainOdf, grainDensityNormalization] = ...
  as_harmonic_odf(grainRbfOdf, harmonicBandwidth);
clear pixelRbfOdf grainRbfOdf

commonGrid = equispacedSO3Grid(pixelOrientations.CS, ...
  pixelOrientations.SS, "resolution", ...
  options.grid_resolution_deg * degree);

[pixelRow, pixelDistribution, pixelDiagnostics] = ...
  one_weighting(pixelOrientations, ...
  ones(length(pixelOrientations),1), pixelIds, "pixel_weighted", ...
  pixelOdf, commonGrid, sampleMeta, options);
[grainRow, grainDistribution, grainDiagnostics] = ...
  one_weighting(grainMeanOrientations, ...
  grainAreas, grainIds, "area_weighted_grain_mean", grainOdf, ...
  commonGrid, sampleMeta, options);
summary = [pixelRow; grainRow];
distribution = [pixelDistribution; grainDistribution];

model = struct();
model.pixel_weighted_odf = pixelOdf;
model.area_weighted_grain_mean_odf = grainOdf;
model.kernel_halfwidth_deg = options.kernel_halfwidth_deg;
model.grid_resolution_deg = options.grid_resolution_deg;
model.harmonic_bandwidth = harmonicBandwidth;
model.common_grid = commonGrid;
model.pixel_weighted_density_normalization_factor = ...
  pixelDensityNormalization;
model.area_weighted_grain_mean_density_normalization_factor = ...
  grainDensityNormalization;
model.pixel_weighted_entropy_diagnostics = pixelDiagnostics;
model.area_weighted_grain_mean_entropy_diagnostics = grainDiagnostics;
end

function [row, distribution, diagnostics] = one_weighting( ...
  orientations, weights, sourceIds, weighting, odf, commonGrid, meta, ...
  options)
cAxis = Miller(0,0,0,1,orientations.CS);
cDirections = orientations * cAxis;
[x, y, z] = canonical_c_axis_components(cDirections, ...
  options.c_axis_component_zero_tolerance);
acuteDeg = acosd(min(1, max(0, x)));
azimuthDeg = mod(atan2d(z,y), 360);
weights = weights(:);
normalizedWeights = weights / sum(weights);

maxMrd = double(max(odf, "resolution", ...
  options.grid_resolution_deg * degree));
textureIndex = double(norm(odf)^2);
mIndex = double(calcMIndex(odf));
[odfEntropy, diagnostics] = nonnegative_grid_entropy(odf, commonGrid);
angleQuantiles = weighted_quantile(acuteDeg, weights, [0.1 0.5 0.9]);
azimuthResultant = abs(sum(normalizedWeights .* exp(1i * deg2rad(azimuthDeg))));

row = table(string(meta.sample), double(meta.diameter_mm), ...
  double(meta.cold_reduction_percent), string(meta.variant), ...
  string(weighting), double(options.kernel_halfwidth_deg), ...
  double(options.grid_resolution_deg), maxMrd, textureIndex, mIndex, ...
  odfEntropy, sum(normalizedWeights .* acuteDeg), angleQuantiles(2), ...
  angleQuantiles(1), angleQuantiles(3), ...
  sum(normalizedWeights(acuteDeg <= 15)), ...
  sum(normalizedWeights(acuteDeg <= 30)), azimuthResultant, ...
  'VariableNames', ["sample", "diameter_mm", ...
  "cold_reduction_percent", "variant", "weighting", ...
  "texture_kernel_halfwidth_deg", "texture_grid_resolution_deg", ...
  "max_mrd", "texture_index", "m_index", "clipped_grid_entropy", ...
  "c_axis_ad_mean_deg", "c_axis_ad_median_deg", ...
  "c_axis_ad_p10_deg", "c_axis_ad_p90_deg", ...
  "c_axis_within_15deg_fraction", "c_axis_within_30deg_fraction", ...
  "c_axis_azimuth_resultant"]);

n = length(orientations);
distribution = table(repmat(string(meta.sample),n,1), ...
  repmat(double(meta.diameter_mm),n,1), ...
  repmat(double(meta.cold_reduction_percent),n,1), ...
  repmat(string(meta.variant),n,1), repmat(string(weighting),n,1), ...
  sourceIds, weights, normalizedWeights, x, y, z, acuteDeg, azimuthDeg, ...
  'VariableNames', ["sample", "diameter_mm", ...
  "cold_reduction_percent", "variant", "weighting", ...
  "source_id", ...
  "source_weight", "normalized_weight", "c_axis_x_ad", ...
  "c_axis_y_td_rd", "c_axis_z_nd", "c_axis_ad_acute_deg", ...
  "c_axis_azimuth_about_ad_deg"]);
end

function [odf, densityNormalization] = ...
  as_harmonic_odf(inputOdf, targetBandwidth)
inputMean = double(mean(inputOdf));
assert(isscalar(inputMean), ...
  "Texture density must be scalar-valued; input ODF size is %s.", ...
  mat2str(size(inputOdf)));
assert(isreal(inputMean) && isfinite(inputMean) && inputMean > 0, ...
  "Texture density has invalid mean %.15g.", inputMean);
densityNormalization = inputMean;
normalizedInputOdf = inputOdf / densityNormalization;
if isa(normalizedInputOdf, "SO3FunHarmonic")
  odf = normalizedInputOdf;
  if odf.bandwidth > targetBandwidth
    odf.bandwidth = targetBandwidth;
  end
else
  odf = SO3FunHarmonic(normalizedInputOdf, ...
    "bandwidth", targetBandwidth);
end
assert(isa(odf, "SO3FunHarmonic"));
assert(odf.bandwidth <= targetBandwidth);
outputMean = double(mean(odf));
assert(isscalar(outputMean) && abs(outputMean - 1) < 1e-6, ...
  ['Harmonic ODF normalization changed during representation conversion: ' ...
  'input_class=%s input_mean=%.15g output_mean=%.15g bandwidth=%g.'], ...
  class(inputOdf), inputMean, outputMean, odf.bandwidth);
end

function [value, diagnostics] = nonnegative_grid_entropy(odf, grid)
rawValues = real(eval(odf, grid)); %#ok<EV2IN>
rawValues = rawValues(:);
assert(~isempty(rawValues) && all(isfinite(rawValues)));
negative = rawValues < 0;
clippedValues = max(rawValues, 0);
normalizationMean = mean(clippedValues);
assert(isfinite(normalizationMean) && normalizationMean > 0);
normalizedDensity = clippedValues / normalizationMean;
positive = normalizedDensity > 0;
value = -sum(normalizedDensity(positive) .* ...
  log(normalizedDensity(positive))) / numel(normalizedDensity);
value = double(value);
assert(isscalar(value) && isreal(value) && isfinite(value));

diagnostics = struct();
diagnostics.grid_point_count = numel(rawValues);
diagnostics.minimum_raw_odf = min(rawValues);
diagnostics.negative_grid_fraction = nnz(negative) / numel(rawValues);
diagnostics.clipped_negative_l1_fraction = ...
  sum(abs(rawValues(negative))) / sum(abs(rawValues));
diagnostics.clipped_normalization_mean = normalizationMean;
end

function ids = source_ids(options, fieldName, n)
if isfield(options, fieldName)
  ids = double(options.(fieldName)(:));
  assert(numel(ids) == n && all(isfinite(ids)) && ...
    numel(unique(ids)) == n, "%s must contain unique IDs.", fieldName);
else
  ids = (1:n)';
end
end

function [x, y, z] = canonical_c_axis_components(directions, tolerance)
x = double(directions.x(:));
y = double(directions.y(:));
z = double(directions.z(:));
lengths = sqrt(x.^2 + y.^2 + z.^2);
assert(all(isfinite(lengths) & lengths > 0));
x = x ./ lengths;
y = y ./ lengths;
z = z ./ lengths;
x(abs(x) < tolerance) = 0;
y(abs(y) < tolerance) = 0;
z(abs(z) < tolerance) = 0;
flip = x < -tolerance | (x == 0 & y < -tolerance) | ...
  (x == 0 & y == 0 & z < 0);
x(flip) = -x(flip);
y(flip) = -y(flip);
z(flip) = -z(flip);
end

function q = weighted_quantile(values, weights, probabilities)
[values, order] = sort(values(:));
weights = weights(order);
cumulative = cumsum(weights) / sum(weights);
q = zeros(size(probabilities));
for k = 1:numel(probabilities)
  index = find(cumulative >= probabilities(k), 1, "first");
  q(k) = values(index);
end
end
