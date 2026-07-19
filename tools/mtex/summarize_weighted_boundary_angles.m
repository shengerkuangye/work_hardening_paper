function [stats, distribution] = summarize_weighted_boundary_angles( ...
  thetaDeg, segLengthUm, binEdgesDeg)
%SUMMARIZE_WEIGHTED_BOUNDARY_ANGLES Summarize boundary angles by segment length.
thetaDeg = double(thetaDeg(:));
segLengthUm = double(segLengthUm(:));
binEdgesDeg = double(binEdgesDeg(:)');
expectedEdges = 1:0.5:94;

assert(~isempty(thetaDeg) && numel(thetaDeg) == numel(segLengthUm), ...
  "thetaDeg and segLengthUm must be nonempty vectors of equal length.");
assert(all(isfinite(thetaDeg)) && all(isfinite(segLengthUm) & segLengthUm > 0), ...
  "Angles must be finite and segment lengths must be positive finite values.");
assert(isequal(binEdgesDeg, expectedEdges), ...
  "binEdgesDeg must equal 1:0.5:94.");
assert(all(thetaDeg >= binEdgesDeg(1) & thetaDeg <= binEdgesDeg(end)), ...
  "Angles fall outside the requested bin range.");

totalLengthUm = sum(segLengthUm);
weightedMean = sum(thetaDeg .* segLengthUm) / totalLengthUm;
[sortedTheta, order] = sort(thetaDeg);
sortedWeight = segLengthUm(order);
weightedMedian = sortedTheta(find(cumsum(sortedWeight) >= ...
  0.5 * totalLengthUm, 1));

[~,~,binIndex] = histcounts(thetaDeg, binEdgesDeg);
binCount = numel(binEdgesDeg) - 1;
boundaryLengthUm = accumarray(binIndex(binIndex > 0), ...
  segLengthUm(binIndex > 0), [binCount,1], @sum, 0);
segmentCount = accumarray(binIndex(binIndex > 0), 1, ...
  [binCount,1], @sum, 0);
binLowerDeg = binEdgesDeg(1:end-1)';
binUpperDeg = binEdgesDeg(2:end)';
binWidthDeg = binUpperDeg - binLowerDeg;
binCenterDeg = (binLowerDeg + binUpperDeg) / 2;
lengthFraction = boundaryLengthUm / totalLengthUm;
probabilityDensityPerDegree = lengthFraction ./ binWidthDeg;
[~, modeIndex] = max(boundaryLengthUm);

stats.segment_count = numel(thetaDeg);
stats.total_boundary_length_um = totalLengthUm;
stats.weighted_mean_angle_deg = weightedMean;
stats.weighted_median_angle_deg = weightedMedian;
stats.min_angle_deg = min(thetaDeg);
stats.max_angle_deg = max(thetaDeg);
stats.mode_bin_center_deg = binCenterDeg(modeIndex);

distribution = table(binLowerDeg, binUpperDeg, binWidthDeg, ...
  binCenterDeg, segmentCount, boundaryLengthUm, lengthFraction, ...
  probabilityDensityPerDegree, cumsum(lengthFraction), 'VariableNames', ...
  {'bin_lower_deg','bin_upper_deg','bin_width_deg','bin_center_deg', ...
  'segment_count','boundary_length_um','length_fraction', ...
  'probability_density_per_degree','cumulative_length_fraction'});
end
