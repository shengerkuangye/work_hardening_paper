function test_comprehensive_maps_morphology_boundaries(scanRoot, outputDir)
%TEST_COMPREHENSIVE_MAPS_MORPHOLOGY_BOUNDARIES Task 3 unit/integration test.

assert(~isempty(which("grain2d")), ...
  "MTEX is not loaded in the current MATLAB session.");
test_synthetic_morphology();
test_synthetic_boundaries();

if nargin >= 1
  assert(nargin == 2, "Full-data testing requires scanRoot and outputDir.");
  test_full_data(string(scanRoot), string(outputDir));
end

fprintf("test_comprehensive_maps_morphology_boundaries passed\n");
end

function test_synthetic_morphology()
grains = synthetic_tiled_grains();
meta = table("synthetic", 1, 0, "raw", 'VariableNames', ...
  {'sample','diameter_mm','cold_reduction_percent','variant'});
options = struct("grain_detection_deg", 2, "min_grain_pixels", 5);
[byGrain, summary, quantiles] = compute_grain_morphology_metrics( ...
  grains, meta, options);

required = ["grain_id","num_pixel","boundary_touching","area_um2", ...
  "ecd_um","perimeter_um","ellipse_long_axis_um", ...
  "ellipse_short_axis_um","aspect_ratio", ...
  "long_axis_ad_angle_deg","max_feret_um","min_feret_um", ...
  "shape_factor"];
assert(all(ismember(required, string(byGrain.Properties.VariableNames))));
assert(isequal(string(summary.Properties.VariableNames), ...
  comprehensive_ebsd_output_contract().summaryColumns. ...
  grain_morphology_summary));

center = byGrain(byGrain.grain_id == 5, :);
assert(height(center) == 1 && ~center.boundary_touching);
assert(abs(center.area_um2 - 2) < 1e-12);
assert(abs(center.ecd_um - sqrt(8 / pi)) < 1e-12);
assert(abs(center.perimeter_um - 6) < 1e-12);
assert(abs(center.max_feret_um - sqrt(5)) < 1e-12);
assert(abs(center.min_feret_um - 1) < 1e-12);
assert(abs(center.shape_factor - 6 / sqrt(8 * pi)) < 1e-12);
assert(abs(center.long_axis_ad_angle_deg) < 1e-10);

[~, ellipseLong, ellipseShort] = fitEllipse(grains(5));
assert(abs(center.ellipse_long_axis_um - 2 * norm(ellipseLong)) < 1e-12);
assert(abs(center.ellipse_short_axis_um - 2 * norm(ellipseShort)) < 1e-12);
assert(abs(center.aspect_ratio - norm(ellipseLong) / ...
  norm(ellipseShort)) < 1e-12);
assert(all(quantiles.long_axis_ad_angle_deg >= 0 & ...
  quantiles.long_axis_ad_angle_deg <= 90 | ...
  isnan(quantiles.long_axis_ad_angle_deg)));
assert(all(ismember(["all","exclude_boundary_touching"], ...
  unique(quantiles.scope))));
assert(all(ismember(["number","area"], unique(quantiles.weighting))));
end

function grains = synthetic_tiled_grains()
% Nine unsmoothed rectangles; grain 5 is a 2 x 1 interior rectangle.
xBreaks = [0 1 3 4];
yBreaks = [0 1 2 3];
[xGrid, yGrid] = meshgrid(xBreaks, yBreaks);
vertices = [xGrid(:), yGrid(:)];
vertexId = reshape(1:numel(xGrid), size(xGrid));
polygons = cell(9, 1);
grainIndex = 0;
for yCell = 1:3
  for xCell = 1:3
    grainIndex = grainIndex + 1;
    ll = vertexId(yCell, xCell);
    lr = vertexId(yCell, xCell + 1);
    ur = vertexId(yCell + 1, xCell + 1);
    ul = vertexId(yCell + 1, xCell);
    polygons{grainIndex} = [ll lr ur ul ll];
  end
end
grains = grain2d(vertices, polygons, [], 'id', (1:9)');
grains.numPixel = repmat(10, 9, 1);
end

function test_synthetic_boundaries()
source = struct();
source.inner_theta_deg = [0.75;1.5;2.5;4;6;14;16];
source.inner_length_um = [1;2;3;4;5;6;7];
source.inner_axis_xyz = repmat([1 0 0], 7, 1);
source.inner_twin_deviation_deg = [80;80;80;80;80;80;0];
source.inner_endpoint_ids = [(1:7)', (2:8)'];
source.outer_theta_deg = [0.75;1.5;2.5;4;6;14;16;85];
source.outer_length_um = [2;3;4;5;6;7;8;9];
source.outer_axis_xyz = repmat([1 0 0], 8, 1);
source.outer_twin_deviation_deg = [80;80;80;80;80;80;20;1];
source.outer_endpoint_ids = [(11:18)', (12:19)'];

meta = table("synthetic", 1, 0, "raw", 'VariableNames', ...
  {'sample','diameter_mm','cold_reduction_percent','variant'});
floors = [0.5 1 2 5];
totalByFloor = nan(size(floors));
for floorIndex = 1:numel(floors)
  options = struct("grain_detection_deg", floors(floorIndex), ...
    "detection_floor_deg", floors(floorIndex), ...
    "classification_deg", 15, "indexed_area_um2", 100, ...
    "min_boundary_axis_deg", 5, ...
    "twin_candidate_tolerance_deg", 5);
  [segments, summary, distribution] = ...
    compute_boundary_network_metrics(source, meta, options);
  assert(isequal(string(summary.Properties.VariableNames), ...
    comprehensive_ebsd_output_contract().summaryColumns.boundary_summary));
  assert(all(segments.misorientation_deg >= floors(floorIndex)));
  assert(all(isnan(segments.axis_x(segments.misorientation_deg < 5))));
  assert(abs(sum(segments.length_um) - ...
    summary.total_boundary_length_um) < 1e-12);
  assert(abs(distribution.cumulative_length_fraction(end) - 1) < 1e-12);
  totalByFloor(floorIndex) = summary.total_boundary_length_um;
end
assert(all(diff(totalByFloor) <= 0));

options.detection_floor_deg = 2;
options.grain_detection_deg = 2;
[segments, summary] = compute_boundary_network_metrics(source, meta, options);
assert(isequal(unique(segments.threshold_bin), ...
  ["2_to_lt5";"5_to_lt15";"ge15"]));
assert(abs(summary.lagb_2_5_length_um + ...
  summary.lagb_5_15_length_um + summary.hagb_ge15_length_um - ...
  summary.total_boundary_length_um) < 1e-12);
assert(summary.lagb_2_5_length_um == 7);
assert(summary.lagb_5_15_length_um == 11);
assert(summary.hagb_ge15_length_um == 17);
assert(summary.total_boundary_length_um == 35);
assert(summary.lagb_2_15_length_fraction == 18 / 35);
assert(summary.hagb_ge15_length_fraction == 17 / 35);
assert(summary.lagb_2_15_count_fraction == 4 / 6);
assert(summary.hagb_ge15_count_fraction == 2 / 6);
assert(summary.twin_candidate_length_fraction == 9 / 35);
assert(all(segments.twin_candidate_label(segments.twin_candidate) == ...
  "angular-axis candidate (non-unique)"));
assert(all(segments.source_class(segments.threshold_bin == "ge15") == ...
  "outer_hagb"));
end

function test_full_data(scanRoot, outputDir)
assert(isfolder(scanRoot), "EBSD scan folder not found: %s", scanRoot);
catalog = comprehensive_ebsd_catalog(scanRoot);
[beforeBytes, beforeTimes] = input_file_stats(catalog.input_path);

[morphology, morphologySummary, boundarySegments, boundarySummary, ...
  boundaryDistribution, sensitivity] = ...
  generate_comprehensive_maps_morphology_boundaries(scanRoot, outputDir);

assert(height(morphology) > 0 && height(boundarySegments) > 0);
assert(height(morphologySummary) == 12);
assert(height(boundarySummary) == 12);
assert(height(sensitivity) == 48);
assert(nnz(morphologySummary.variant == "raw") == 6);
assert(nnz(morphologySummary.variant == "denoised") == 6);
assert(nnz(boundarySummary.variant == "raw") == 6);
assert(nnz(boundarySummary.variant == "denoised") == 6);
assert(isequal(unique(sensitivity.grain_detection_deg), ...
  [0.5;1;2;5]));
floorGroups = findgroups(sensitivity.sample, sensitivity.variant);
assert(all(splitapply(@numel, sensitivity.grain_detection_deg, ...
  floorGroups) == 4));
assert(all(boundarySegments.grain_detection_deg == 2));
assert(all(boundarySegments.detection_floor_deg == 2));
assert(all(boundaryDistribution.grain_detection_deg == 2));
assert(all(boundaryDistribution.detection_floor_deg == 2));
assert(all(boundarySummary.total_boundary_length_um > 0));
assert(all(abs(boundarySummary.lagb_2_5_length_um + ...
  boundarySummary.lagb_5_15_length_um + ...
  boundarySummary.hagb_ge15_length_um - ...
  boundarySummary.total_boundary_length_um) < 1e-8));
assert(all(boundaryDistribution.cumulative_length_fraction > 0 & ...
  boundaryDistribution.cumulative_length_fraction <= 1));

requiredFiles = [ ...
  fullfile(outputDir, "02_grain_morphology", ...
    "grain_morphology_by_grain.csv")
  fullfile(outputDir, "02_grain_morphology", ...
    "grain_morphology_summary.csv")
  fullfile(outputDir, "02_grain_morphology", ...
    "grain_morphology_quantiles.csv")
  fullfile(outputDir, "02_grain_morphology", ...
    "grain_morphology_trends.png")
  fullfile(outputDir, "03_boundaries", "boundary_segments.csv")
  fullfile(outputDir, "03_boundaries", "boundary_summary.csv")
  fullfile(outputDir, "03_boundaries", ...
    "boundary_angle_distributions.csv")
  fullfile(outputDir, "03_boundaries", ...
    "boundary_detection_sensitivity.csv")
  fullfile(outputDir, "03_boundaries", "boundary_trends.png")
];
for filePath = requiredFiles'
  assert(isfile(filePath), "Missing Task 3 artifact: %s", filePath);
end
for catalogIndex = 1:height(catalog)
  mapFile = fullfile(outputDir, "01_standard_maps", ...
    catalog.sample(catalogIndex) + "_" + ...
    catalog.variant(catalogIndex) + "_maps.png");
  assert(isfile(mapFile), "Missing standard-map panel: %s", mapFile);
  mapInfo = dir(mapFile);
  assert(mapInfo.bytes > 0, "Empty standard-map panel: %s", mapFile);
end

morphologyCsv = readtable(fullfile(outputDir, "02_grain_morphology", ...
  "grain_morphology_summary.csv"));
boundaryCsv = readtable(fullfile(outputDir, "03_boundaries", ...
  "boundary_summary.csv"));
sensitivityCsv = readtable(fullfile(outputDir, "03_boundaries", ...
  "boundary_detection_sensitivity.csv"));
assert(isequal(string(morphologyCsv.Properties.VariableNames), ...
  comprehensive_ebsd_output_contract().summaryColumns. ...
  grain_morphology_summary));
assert(isequal(string(boundaryCsv.Properties.VariableNames), ...
  comprehensive_ebsd_output_contract().summaryColumns.boundary_summary));
assert(isequal(string(sensitivityCsv.Properties.VariableNames), ...
  comprehensive_ebsd_output_contract().summaryColumns.boundary_summary));
assert(height(morphologyCsv) == 12 && height(boundaryCsv) == 12 && ...
  height(sensitivityCsv) == 48);

[afterBytes, afterTimes] = input_file_stats(catalog.input_path);
assert(isequal(afterBytes, beforeBytes));
assert(isequal(afterTimes, beforeTimes));
end

function [bytes, times] = input_file_stats(paths)
bytes = nan(numel(paths), 1);
times = nan(numel(paths), 1);
for pathIndex = 1:numel(paths)
  info = dir(paths(pathIndex));
  assert(isscalar(info));
  bytes(pathIndex) = info.bytes;
  times(pathIndex) = info.datenum;
end
end
