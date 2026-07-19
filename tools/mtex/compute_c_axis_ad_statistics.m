function [stats, distribution] = compute_c_axis_ad_statistics( ...
  thetaDeg, totalPixelCount, binEdgesDeg)
%COMPUTE_C_AXIS_AD_STATISTICS Summarize acute c-axis angles to bar AD.

thetaDeg = double(thetaDeg(:));
binEdgesDeg = double(binEdgesDeg(:)');
n = numel(thetaDeg);

assert(n > 0, "thetaDeg must not be empty.");
assert(all(isfinite(thetaDeg)), "thetaDeg contains nonfinite values.");
assert(all(thetaDeg >= 0 & thetaDeg <= 90), ...
  "thetaDeg values must be between 0 and 90 degrees.");
assert(isscalar(totalPixelCount) && isfinite(totalPixelCount) && ...
  totalPixelCount > 0 && totalPixelCount >= n && ...
  totalPixelCount == fix(totalPixelCount), ...
  "totalPixelCount must be a positive integer not less than the valid pixel count.");
assert(numel(binEdgesDeg) >= 2 && binEdgesDeg(1) == 0 && ...
  binEdgesDeg(end) == 90 && all(diff(binEdgesDeg) > 0), ...
  "binEdgesDeg must increase from 0 to 90 degrees.");

cos2 = cosd(thetaDeg).^2;
stats.total_pixel_count = totalPixelCount;
stats.alpha_ti_pixel_count = n;
stats.alpha_ti_pixel_fraction = n / totalPixelCount;
stats.mean_angle_deg = mean(thetaDeg);
stats.std_angle_deg = std(thetaDeg, 1);
quartiles = prctile(thetaDeg, [25, 50, 75]);
stats.q25_angle_deg = quartiles(1);
stats.median_angle_deg = quartiles(2);
stats.q75_angle_deg = quartiles(3);
stats.mean_cos2 = mean(cos2);
stats.alignment_factor = (3 * stats.mean_cos2 - 1) / 2;
stats.fraction_within_15deg = mean(thetaDeg <= 15);
stats.fraction_within_30deg = mean(thetaDeg <= 30);
stats.fraction_within_45deg = mean(thetaDeg <= 45);

pixelCount = histcounts(thetaDeg, binEdgesDeg)';
binLowerDeg = binEdgesDeg(1:end-1)';
binUpperDeg = binEdgesDeg(2:end)';
binWidthDeg = binUpperDeg - binLowerDeg;
binCenterDeg = (binLowerDeg + binUpperDeg) / 2;
probability = pixelCount / n;
probabilityDensityPerDegree = probability ./ binWidthDeg;

distribution = table(binLowerDeg, binUpperDeg, binWidthDeg, ...
  binCenterDeg, pixelCount, probability, probabilityDensityPerDegree, ...
  'VariableNames', {'bin_lower_deg', 'bin_upper_deg', 'bin_width_deg', ...
  'bin_center_deg', 'pixel_count', 'probability', ...
  'probability_density_per_degree'});
end
