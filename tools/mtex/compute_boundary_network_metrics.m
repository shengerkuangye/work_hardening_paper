function [segments, summary, distribution, audit] = ...
  compute_boundary_network_metrics(source, meta, options)
%COMPUTE_BOUNDARY_NETWORK_METRICS Length-weighted unsmoothed Ti-Hex network.

arguments
  source (1,1) struct
  meta table
  options (1,1) struct
end
requiredSource = ["inner_theta_deg","inner_length_um", ...
  "outer_theta_deg","outer_length_um"];
assert(all(isfield(source, requiredSource)));
assert(height(meta) == 1 && all(ismember( ...
  ["sample","diameter_mm","cold_reduction_percent","variant"], ...
  string(meta.Properties.VariableNames))));
requiredOptions = ["grain_detection_deg","detection_floor_deg", ...
  "classification_deg","indexed_area_um2","min_boundary_axis_deg", ...
  "twin_candidate_tolerance_deg"];
assert(all(isfield(options, requiredOptions)));

floorDeg = double(options.detection_floor_deg);
classificationDeg = double(options.classification_deg);
assert(isscalar(floorDeg) && isfinite(floorDeg) && floorDeg >= 0);
assert(isscalar(classificationDeg) && isfinite(classificationDeg) && ...
  classificationDeg > max(floorDeg, 5));

inner = normalize_source_part(source, "inner");
outer = normalize_source_part(source, "outer");
inner = select_source_rows(inner, inner.theta_deg >= floorDeg);
outer = select_source_rows(outer, outer.theta_deg >= floorDeg);

[~, ~, audit, masks] = partition_ti_hex_boundary_segments( ...
  inner.theta_deg, inner.length_um, outer.theta_deg, outer.length_um, ...
  floorDeg, classificationDeg);
innerEligible = select_source_rows(inner, masks.inner_lagb);
outerEligible = select_source_rows(outer, masks.outer_hagb);
eligible = concatenate_sources(innerEligible, outerEligible);
sourceClass = [repmat("inner_lagb", numel(innerEligible.theta_deg), 1); ...
  repmat("outer_hagb", numel(outerEligible.theta_deg), 1)];
assert(numel(eligible.theta_deg) == audit.eligible_segment_count);
assert(abs(sum(eligible.length_um) - ...
  audit.total_eligible_boundary_length_um) < 1e-8);

stableAxis = eligible.theta_deg >= options.min_boundary_axis_deg & ...
  all(isfinite(eligible.axis_xyz), 2);
eligible.axis_xyz(~stableAxis, :) = NaN;
twinCandidate = stableAxis & ...
  isfinite(eligible.twin_deviation_deg) & ...
  eligible.twin_deviation_deg <= options.twin_candidate_tolerance_deg;
twinLabel = repmat("not candidate", numel(twinCandidate), 1);
twinLabel(twinCandidate) = "angular-axis candidate (non-unique)";
thresholdBin = classify_threshold_bins(eligible.theta_deg);

segment_id = (1:numel(eligible.theta_deg))';
misorientation_deg = eligible.theta_deg;
length_um = eligible.length_um;
source_class = sourceClass;
threshold_bin = thresholdBin;
axis_stable = stableAxis;
axis_x = eligible.axis_xyz(:, 1);
axis_y = eligible.axis_xyz(:, 2);
axis_z = eligible.axis_xyz(:, 3);
twin_candidate_deviation_deg = eligible.twin_deviation_deg;
twin_candidate = twinCandidate;
twin_candidate_label = twinLabel;
endpoint_ebsd_id_1 = eligible.endpoint_ids(:, 1);
endpoint_ebsd_id_2 = eligible.endpoint_ids(:, 2);
segments = table(segment_id, misorientation_deg, length_um, ...
  source_class, threshold_bin, axis_stable, axis_x, axis_y, axis_z, ...
  twin_candidate_deviation_deg, twin_candidate, twin_candidate_label, ...
  endpoint_ebsd_id_1, endpoint_ebsd_id_2);
segments = add_common_columns(segments, meta);
segments.grain_detection_deg = repmat( ...
  double(options.grain_detection_deg), height(segments), 1);
segments.detection_floor_deg = repmat(floorDeg, height(segments), 1);
segments.classification_deg = repmat(classificationDeg, height(segments), 1);

stats = calculate_boundary_threshold_metrics(eligible.theta_deg, ...
  eligible.length_um, floorDeg, classificationDeg, ...
  options.indexed_area_um2);
lagb2to5 = thresholdBin == "2_to_lt5";
lagb5to15 = thresholdBin == "5_to_lt15";
hagb = thresholdBin == "ge15";
lagb2to15 = lagb2to5 | lagb5to15;
totalLength = stats.total_eligible_boundary_length_um;
totalCount = stats.eligible_segment_count;
twinLength = sum(eligible.length_um(twinCandidate));

summary = table(string(meta.sample), double(meta.diameter_mm), ...
  double(meta.cold_reduction_percent), string(meta.variant), ...
  double(options.grain_detection_deg), classificationDeg, totalLength, ...
  stats.total_eligible_length_density_um_per_um2, ...
  sum(eligible.length_um(lagb2to5)), ...
  sum(eligible.length_um(lagb5to15)), ...
  sum(eligible.length_um(hagb)), ...
  sum(eligible.length_um(lagb2to15)) / totalLength, ...
  sum(eligible.length_um(hagb)) / totalLength, ...
  nnz(lagb2to15) / totalCount, nnz(hagb) / totalCount, ...
  twinLength / totalLength, 'VariableNames', cellstr( ...
  comprehensive_ebsd_output_contract().summaryColumns.boundary_summary));

distribution = cumulative_distribution(eligible.theta_deg, ...
  eligible.length_um, meta, options);
audit.bin_length_sum_um = sum(eligible.length_um( ...
  thresholdBin == "below2" | thresholdBin == "2_to_lt5" | ...
  thresholdBin == "5_to_lt15" | thresholdBin == "ge15"));
assert(abs(audit.bin_length_sum_um - totalLength) < 1e-8);
end

function part = normalize_source_part(source, prefix)
part.theta_deg = double(source.(prefix + "_theta_deg")(:));
part.length_um = double(source.(prefix + "_length_um")(:));
n = numel(part.theta_deg);
assert(numel(part.length_um) == n && ...
  all(isfinite(part.theta_deg)) && ...
  all(isfinite(part.length_um) & part.length_um > 0));
axisName = prefix + "_axis_xyz";
deviationName = prefix + "_twin_deviation_deg";
endpointName = prefix + "_endpoint_ids";
if isfield(source, axisName)
  part.axis_xyz = double(source.(axisName));
else
  part.axis_xyz = nan(n, 3);
end
if isfield(source, deviationName)
  part.twin_deviation_deg = double(source.(deviationName)(:));
else
  part.twin_deviation_deg = nan(n, 1);
end
if isfield(source, endpointName)
  part.endpoint_ids = double(source.(endpointName));
else
  part.endpoint_ids = nan(n, 2);
end
assert(isequal(size(part.axis_xyz), [n 3]));
assert(numel(part.twin_deviation_deg) == n);
assert(isequal(size(part.endpoint_ids), [n 2]));
end

function output = select_source_rows(input, mask)
fields = string(fieldnames(input));
for field = fields'
  output.(field) = input.(field)(mask, :);
end
end

function output = concatenate_sources(inner, outer)
fields = string(fieldnames(inner));
for field = fields'
  output.(field) = [inner.(field); outer.(field)];
end
end

function bins = classify_threshold_bins(thetaDeg)
bins = strings(size(thetaDeg));
bins(thetaDeg < 2) = "below2";
bins(thetaDeg >= 2 & thetaDeg < 5) = "2_to_lt5";
bins(thetaDeg >= 5 & thetaDeg < 15) = "5_to_lt15";
bins(thetaDeg >= 15) = "ge15";
assert(all(strlength(bins) > 0));
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

function distribution = cumulative_distribution(thetaDeg, lengthUm, ...
  meta, options)
[angleValues, ~, group] = unique(thetaDeg);
segmentCount = accumarray(group, 1);
boundaryLengthUm = accumarray(group, lengthUm);
cumulativeLengthFraction = cumsum(boundaryLengthUm) / sum(lengthUm);
cumulativeNumberFraction = cumsum(segmentCount) / numel(thetaDeg);
n = numel(angleValues);
distribution = table(repmat(string(meta.sample), n, 1), ...
  repmat(double(meta.diameter_mm), n, 1), ...
  repmat(double(meta.cold_reduction_percent), n, 1), ...
  repmat(string(meta.variant), n, 1), ...
  repmat(double(options.grain_detection_deg), n, 1), ...
  repmat(double(options.detection_floor_deg), n, 1), ...
  angleValues, segmentCount, boundaryLengthUm, ...
  cumulativeLengthFraction, cumulativeNumberFraction, ...
  'VariableNames', {'sample','diameter_mm','cold_reduction_percent', ...
  'variant','grain_detection_deg','detection_floor_deg', ...
  'misorientation_deg','segment_count','boundary_length_um', ...
  'cumulative_length_fraction','cumulative_number_fraction'});
end
