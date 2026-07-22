function colorLimits = comprehensive_ebsd_change_map_limits( ...
  metricMinimum, metricMaximum)
%COMPREHENSIVE_EBSD_CHANGE_MAP_LIMITS Set shared limits by map metric.

assert(isnumeric(metricMinimum) && isnumeric(metricMaximum));
assert(isequal(size(metricMinimum), size(metricMaximum)) && ...
  size(metricMinimum, 2) == 6, ...
  "Metric extrema must be equally sized arrays with six columns.");

colorLimits = nan(6, 2);
colorLimits(1:2, :) = repmat([0 1], 2, 1);
orientationMaximum = max(metricMaximum(:, 3), [], "omitnan");
if isempty(orientationMaximum) || ~isfinite(orientationMaximum) || ...
    orientationMaximum <= 0
  orientationMaximum = 1;
end
colorLimits(3, :) = [0 orientationMaximum];
for metricIndex = 4:6
  extrema = [metricMinimum(:, metricIndex); ...
    metricMaximum(:, metricIndex)];
  scale = max(abs(extrema), [], "omitnan");
  if isempty(scale) || ~isfinite(scale) || scale <= 0
    scale = 1;
  end
  colorLimits(metricIndex, :) = [-scale scale];
end
end
