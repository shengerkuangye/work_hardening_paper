function stats = calculate_boundary_threshold_metrics(thetaDeg, ...
  segLengthUm, detectionFloorDeg, classificationAngleDeg, indexedAreaUm2)
%CALCULATE_BOUNDARY_THRESHOLD_METRICS Length-weighted LAGB/HAGB metrics.

thetaDeg = double(thetaDeg(:));
segLengthUm = double(segLengthUm(:));

assert(~isempty(thetaDeg) && numel(thetaDeg) == numel(segLengthUm), ...
  "Angles and segment lengths must be nonempty vectors of equal length.");
assert(all(isfinite(thetaDeg)) && ...
  all(isfinite(segLengthUm) & segLengthUm > 0), ...
  "Angles must be finite and segment lengths must be positive finite values.");
assert(isscalar(detectionFloorDeg) && isfinite(detectionFloorDeg) && ...
  detectionFloorDeg >= 0, ...
  "The detection floor must be a nonnegative finite scalar.");
assert(isscalar(classificationAngleDeg) && ...
  isfinite(classificationAngleDeg) && ...
  classificationAngleDeg > detectionFloorDeg, ...
  "The classification angle must exceed the detection floor.");
assert(all(thetaDeg >= detectionFloorDeg), ...
  "Angles below the detection floor are not eligible.");
assert(isscalar(indexedAreaUm2) && isfinite(indexedAreaUm2) && ...
  indexedAreaUm2 > 0, "Indexed area must be a positive finite scalar.");

lagbMask = thetaDeg < classificationAngleDeg;
hagbMask = ~lagbMask;
totalLengthUm = sum(segLengthUm);
lagbLengthUm = sum(segLengthUm(lagbMask));
hagbLengthUm = sum(segLengthUm(hagbMask));

[sortedThetaDeg, order] = sort(thetaDeg);
sortedLengthUm = segLengthUm(order);
weightedMedianIndex = find(cumsum(sortedLengthUm) >= ...
  0.5 * totalLengthUm, 1);

stats.eligible_segment_count = numel(thetaDeg);
stats.total_eligible_boundary_length_um = totalLengthUm;
stats.lagb_segment_count = nnz(lagbMask);
stats.hagb_segment_count = nnz(hagbMask);
stats.lagb_length_um = lagbLengthUm;
stats.hagb_length_um = hagbLengthUm;
stats.lagb_length_fraction = lagbLengthUm / totalLengthUm;
stats.hagb_length_fraction = hagbLengthUm / totalLengthUm;
stats.total_eligible_length_density_um_per_um2 = ...
  totalLengthUm / indexedAreaUm2;
stats.lagb_length_density_um_per_um2 = lagbLengthUm / indexedAreaUm2;
stats.hagb_length_density_um_per_um2 = hagbLengthUm / indexedAreaUm2;
stats.weighted_mean_angle_deg = ...
  sum(thetaDeg .* segLengthUm) / totalLengthUm;
stats.weighted_median_angle_deg = sortedThetaDeg(weightedMedianIndex);
stats.min_angle_deg = min(thetaDeg);
stats.max_angle_deg = max(thetaDeg);

assert(abs(lagbLengthUm + hagbLengthUm - totalLengthUm) < 1e-8);
assert(abs(stats.lagb_length_fraction + stats.hagb_length_fraction - 1) ...
  < 1e-12);
end
