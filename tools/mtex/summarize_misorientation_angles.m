function [stats, distribution] = summarize_misorientation_angles( ...
  thetaDeg, binEdgesDeg)
%SUMMARIZE_MISORIENTATION_ANGLES Summarize angles using fixed bins.
thetaDeg = double(thetaDeg(:));
binEdgesDeg = double(binEdgesDeg(:)');
assert(~isempty(thetaDeg) && all(isfinite(thetaDeg)), ...
  "thetaDeg must contain finite values.");
assert(all(thetaDeg >= binEdgesDeg(1) & thetaDeg <= binEdgesDeg(end)), ...
  "Angles fall outside the requested bin range.");
assert(binEdgesDeg(1) == 0 && binEdgesDeg(end) == 94 && ...
  all(abs(diff(binEdgesDeg) - 0.1) < 1e-12), ...
  "binEdgesDeg must span 0 to 94 degrees in 0.1-degree bins.");

pairCount = histcounts(thetaDeg, binEdgesDeg)';
binLowerDeg = binEdgesDeg(1:end-1)';
binUpperDeg = binEdgesDeg(2:end)';
binWidthDeg = binUpperDeg - binLowerDeg;
binCenterDeg = (binLowerDeg + binUpperDeg) / 2;
probability = pairCount / numel(thetaDeg);
probabilityDensityPerDegree = probability ./ binWidthDeg;
cumulativeProbability = cumsum(probability);
[~, modeIndex] = max(pairCount);
percentiles = prctile(thetaDeg, [10,25,50,75,90]);

stats.pair_count = numel(thetaDeg);
stats.mean_angle_deg = mean(thetaDeg);
stats.std_angle_deg = std(thetaDeg,1);
stats.q10_angle_deg = percentiles(1);
stats.q25_angle_deg = percentiles(2);
stats.median_angle_deg = percentiles(3);
stats.q75_angle_deg = percentiles(4);
stats.q90_angle_deg = percentiles(5);
stats.mode_bin_center_deg = binCenterDeg(modeIndex);
stats.min_angle_deg = min(thetaDeg);
stats.max_angle_deg = max(thetaDeg);

distribution = table(binLowerDeg, binUpperDeg, binWidthDeg, ...
  binCenterDeg, pairCount, probability, probabilityDensityPerDegree, ...
  cumulativeProbability, 'VariableNames', ...
  {'bin_lower_deg','bin_upper_deg','bin_width_deg','bin_center_deg', ...
  'pair_count','probability','probability_density_per_degree', ...
  'cumulative_probability'});
end
