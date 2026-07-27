function [normalizedDensity,normalizationFactor] = ...
  normalize_positive_mean_density(inputDensity)
%NORMALIZE_POSITIVE_MEAN_DENSITY Normalize scalar texture density to MRD.

inputMean = double(mean(inputDensity));
assert(isscalar(inputMean), ...
  "Texture density mean must be scalar.");
assert(isreal(inputMean) && isfinite(inputMean) && inputMean > 0, ...
  "Texture density mean must be real, finite and positive.");

normalizationFactor = inputMean;
normalizedDensity = inputDensity / normalizationFactor;
outputMean = double(mean(normalizedDensity));
assert(isscalar(outputMean) && isreal(outputMean) && ...
  isfinite(outputMean) && abs(outputMean - 1) < 1e-6, ...
  "Texture density failed to normalize to unit mean.");
end
