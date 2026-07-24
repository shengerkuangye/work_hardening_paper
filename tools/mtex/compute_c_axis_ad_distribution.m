function [distribution,audit] = compute_c_axis_ad_distribution( ...
  thetaDeg,weights,binEdgesDeg)
%COMPUTE_C_AXIS_AD_DISTRIBUTION Axial PDF and random-normalized MRD.

thetaDeg = double(thetaDeg(:));
weights = double(weights(:));
binEdgesDeg = double(binEdgesDeg(:)');
assert(~isempty(thetaDeg) && numel(thetaDeg) == numel(weights));
assert(all(isfinite(thetaDeg) & thetaDeg >= 0 & thetaDeg <= 90));
assert(all(isfinite(weights) & weights > 0));
assert(numel(binEdgesDeg) >= 2 && binEdgesDeg(1) == 0 && ...
  binEdgesDeg(end) == 90 && all(diff(binEdgesDeg) > 0));

binCount = numel(binEdgesDeg)-1;
binIndex = discretize(thetaDeg,binEdgesDeg);
assert(all(~isnan(binIndex)));
observedWeight = accumarray(binIndex,weights,[binCount,1],@sum,0);
totalWeight = sum(weights);
observedProbability = observedWeight / totalWeight;
binLowerDeg = binEdgesDeg(1:end-1)';
binUpperDeg = binEdgesDeg(2:end)';
binWidthDeg = binUpperDeg-binLowerDeg;
binCenterDeg = (binLowerDeg+binUpperDeg)/2;
pdfPerDegree = observedProbability ./ binWidthDeg;
randomProbability = cosd(binLowerDeg)-cosd(binUpperDeg);
assert(all(randomProbability > 0));
mrd = observedProbability ./ randomProbability;

distribution = table(binLowerDeg,binUpperDeg,binCenterDeg, ...
  binWidthDeg,observedWeight,observedProbability,pdfPerDegree, ...
  randomProbability,mrd, ...
  'VariableNames',["bin_lower_deg","bin_upper_deg", ...
  "bin_center_deg","bin_width_deg","observed_weight", ...
  "observed_probability","pdf_per_degree", ...
  "random_probability","mrd"]);
audit = struct();
audit.valid_source_count = numel(thetaDeg);
audit.valid_source_weight = totalWeight;
audit.probability_sum = sum(observedProbability);
audit.random_probability_sum = sum(randomProbability);
assert(abs(audit.probability_sum-1) < 1e-12);
assert(abs(audit.random_probability_sum-1) < 1e-12);
end
