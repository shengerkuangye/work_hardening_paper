function test_comprehensive_intragranular_texture(scanRoot, outputDir)
%TEST_COMPREHENSIVE_INTRAGRANULAR_TEXTURE Verify Task 4 calculations.

arguments
  scanRoot (1,1) string = ""
  outputDir (1,1) string = ""
end

assert(~isempty(which("EBSD")), ...
  "MTEX is not loaded in the current MATLAB session.");
contract = comprehensive_ebsd_output_contract();
assert(contract.parameters.coordinate_x == "AD");
assert(contract.parameters.coordinate_y == "TD_RD");
assert(contract.parameters.coordinate_z == "ND");

cs = crystalSymmetry("6/mmm", [2.95 2.95 4.68], ...
  "mineral", "Ti-Hex");
rotations = repmat(rotation.id, 5, 5);
ebsd = EBSDsquare([], rotations, 2 * ones(5), 1:2, ...
  {"not indexed", cs}, "dxy", [1 1]);
[grains, ebsd.grainId] = calcGrains(ebsd, ...
  "threshold", contract.parameters.primary_grain_detection_deg * degree);

meta = struct("sample", "synthetic", "diameter_mm", 1, ...
  "cold_reduction_percent", 0, "variant", "raw");
intraOptions = struct( ...
  "grain_detection_deg", contract.parameters.primary_grain_detection_deg, ...
  "min_grain_pixels", 1, ...
  "kam_orders", contract.parameters.kam_orders, ...
  "kam_thresholds_deg", contract.parameters.kam_thresholds_deg, ...
  "axis_min_grod_deg", 0.1);
[intraSummary, pixelTable, grainTable] = ...
  compute_intragranular_metrics(ebsd, grains, meta, intraOptions);

assert(isequal(string(intraSummary.Properties.VariableNames), ...
  contract.summaryColumns.intragranular_summary));
assert(height(intraSummary) == 4);
assert(all(intraSummary.valid_pixel_count == 25));
assert(all(abs(intraSummary.kam_mean_deg) < 1e-12));
assert(all(abs(intraSummary.kam_median_deg) < 1e-12));
assert(all(abs(intraSummary.kam_p90_deg) < 1e-12));
assert(all(abs(intraSummary.grod_mean_deg) < 1e-12));
assert(all(abs(intraSummary.grod_median_deg) < 1e-12));
assert(all(abs(intraSummary.grod_p90_deg) < 1e-12));
assert(all(abs(intraSummary.gos_number_mean_deg) < 1e-12));
assert(all(abs(intraSummary.gos_area_weighted_mean_deg) < 1e-12));
assert(height(pixelTable) == 25);
assert(all(pixelTable.grod_angle_deg == 0));
assert(all(pixelTable.grod_axis_valid == 0));
assert(height(grainTable) == 1 && grainTable.gos_deg == 0);

% Real scans contain a sparse cubic phase. The Ti-Hex metrics must remain
% valid without requesting the single-phase-only EBSD.orientations property.
cubicCs = crystalSymmetry("m-3m", [3.3 3.3 3.3], ...
  "mineral", "Titanium cubic");
multiPhaseId = 2 * ones(3);
multiPhaseId(2,2) = 3;
multiEbsd = EBSDsquare([], repmat(rotation.id, 3, 3), multiPhaseId, ...
  1:3, {"not indexed", cs, cubicCs}, "dxy", [1 1]);
[multiGrains, multiEbsd.grainId] = calcGrains(multiEbsd, ...
  "threshold", contract.parameters.primary_grain_detection_deg * degree);
[multiSummary, multiPixel, ~] = compute_intragranular_metrics( ...
  multiEbsd, multiGrains, meta, intraOptions);
assert(height(multiPixel) == 8);
assert(all(isfinite(multiSummary{:, 9:end}), "all"));

% Independent KAM oracle: a 3x3 column gradient [0,3,6] degrees gives the
% center four-neighbor differences [3,3,0,0] and eight-neighbor differences
% [3,3,3,3,3,3,0,0]. These are hand-computed expectations, not another KAM
% call. The 2 degree cutoff retains only the two zero-difference neighbors.
columnAngles = repmat([0 3 6] * degree, 3, 1);
columnRotations = rotation.byAxisAngle(zvector, columnAngles);
columnGradient = EBSDsquare([], columnRotations, 2 * ones(3), 1:2, ...
  {"not indexed", cs}, "dxy", [1 1]);
[columnGrains, columnGradient.grainId] = calcGrains(columnGradient, ...
  "threshold", 5 * degree);
assert(length(columnGrains) == 1);
center = columnGradient.x == 1 & columnGradient.y == 1;
assert(nnz(center) == 1);
dx = abs(double(columnGradient.x - columnGradient.x(center)));
dy = abs(double(columnGradient.y - columnGradient.y(center)));
order1Neighbors = (dx + dy == 1);
order2Neighbors = max(dx,dy) == 1;
assert(nnz(order1Neighbors) == 4 && nnz(order2Neighbors) == 8);
assert(all(columnGradient.grainId(order2Neighbors) == ...
  columnGradient.grainId(center)));
kamOracleOptions = intraOptions;
kamOracleOptions.kam_orders = [1 2];
kamOracleOptions.kam_thresholds_deg = [2 5];
[~, columnPixel, ~] = compute_intragranular_metrics( ...
  columnGradient, columnGrains, meta, kamOracleOptions);
centerRow = columnPixel.x_um == 1 & columnPixel.y_um == 1;
assert(nnz(centerRow) == 1);
kamNumericalToleranceDeg = 2e-6;
assert(abs(columnPixel.kam_order1_threshold5deg(centerRow) - 1.5) < ...
  kamNumericalToleranceDeg);
assert(abs(columnPixel.kam_order2_threshold5deg(centerRow) - 2.25) < ...
  kamNumericalToleranceDeg);
assert(abs(columnPixel.kam_order1_threshold2deg(centerRow)) < ...
  kamNumericalToleranceDeg);
assert(abs(columnPixel.kam_order2_threshold2deg(centerRow)) < ...
  kamNumericalToleranceDeg);

% gridify assigns new regular-grid IDs and preserves source IDs in oldId.
% Deliberately permute a nonuniform grid so matching source IDs to new IDs
% is demonstrably wrong; the exported per-pixel KAM must follow oldId.
gradientAngles = reshape((0:8) * 0.4 * degree, 3, 3);
gradientRotations = rotation.byAxisAngle(zvector, gradientAngles);
orderedGradient = EBSDsquare([], gradientRotations, 2 * ones(3), 1:2, ...
  {"not indexed", cs}, "dxy", [1 1]);
permutation = [9 1 5 3 7 2 8 4 6];
permutedGradient = EBSD(orderedGradient);
permutedGradient.id = permutation(:);
[permutedGrains, permutedGradient.grainId] = calcGrains( ...
  permutedGradient, "threshold", 5 * degree);
mappingOptions = intraOptions;
mappingOptions.kam_orders = 1;
mappingOptions.kam_thresholds_deg = 5;
[~, mappedPixel, mappedGrain] = compute_intragranular_metrics( ...
  permutedGradient, permutedGrains, meta, mappingOptions);
oracleGrid = permutedGradient.gridify;
oracleKam = double(oracleGrid.KAM("order", 1, ...
  "threshold", 5 * degree) / degree);
[foundOldId, oraclePosition] = ismember(mappedPixel.ebsd_id, ...
  double(oracleGrid.oldId(:)));
[~, wrongNewIdPosition] = ismember(mappedPixel.ebsd_id, ...
  double(oracleGrid.id(:)));
assert(all(foundOldId) && any(oraclePosition ~= wrongNewIdPosition));
expectedKam = oracleKam(oraclePosition);
actualKam = mappedPixel.kam_order1_threshold5deg;
assert(isequal(isnan(actualKam), isnan(expectedKam)));
finiteKam = isfinite(expectedKam);
assert(max(abs(actualKam(finiteKam) - expectedKam(finiteKam))) < 1e-12, ...
  "Per-pixel KAM is not mapped through gridify oldId.");
assert(any(actualKam(finiteKam) > 0));
expectedGrod = double(angle(permutedGradient.orientations, ...
  permutedGrains.meanOrientation) / degree);
nominalGrod = [1.6; 1.2; 0.8; 0.4; 0; 0.4; 0.8; 1.2; 1.6];
assert(nnz(expectedGrod > 0.1) == 8);
assert(max(abs(expectedGrod - nominalGrod)) < 1e-5);
assert(max(abs(mappedPixel.grod_angle_deg - expectedGrod)) < 1e-12);
assert(height(mappedGrain) == 1);
expectedGos = double(permutedGrains.GOS / degree);
assert(expectedGos > 0.1);
assert(abs(expectedGos - mean(nominalGrod)) < 1e-5);
assert(abs(mappedGrain.gos_deg - expectedGos) < 1e-12);
assert(mappedGrain.valid_grod_axis_pixel_count == 8);
assert(abs(abs(mappedGrain.grod_crystal_axis_principal_z) - 1) < 1e-8);
assert(abs(abs(mappedGrain.grod_specimen_axis_principal_z_nd) - 1) < 1e-8);

cAlongAd = orientation.byAxisAngle(yvector, 90 * degree, cs);
cAlongMinusAd = orientation.byAxisAngle(yvector, -90 * degree, cs);
cAlongTd = orientation.byAxisAngle(xvector, -90 * degree, cs);
cAlongMinusTd = orientation.byAxisAngle(xvector, 90 * degree, cs);
cAlongNd = orientation.id(cs);
cAlongMinusNd = orientation.byAxisAngle(xvector, pi, cs);
canonicalOrientations = [cAlongAd; cAlongMinusAd; cAlongTd; ...
  cAlongMinusTd; cAlongNd; cAlongMinusNd];
textureOptions = struct( ...
  "kernel_halfwidth_deg", ...
    contract.parameters.texture_kernel_halfwidth_deg, ...
  "grid_resolution_deg", ...
    contract.parameters.texture_grid_resolution_deg, ...
  "c_axis_component_zero_tolerance", ...
    contract.parameters.c_axis_component_zero_tolerance, ...
  "min_grain_pixels", contract.parameters.min_grain_pixels, ...
  "pixel_ids", (1:6)', "grain_ids", (1:6)');
[textureSummary, cAxisDistribution, textureModel] = ...
  compute_texture_metrics(canonicalOrientations, canonicalOrientations, ...
  ones(6,1), meta, textureOptions);

assert(isequal(string(textureSummary.Properties.VariableNames), ...
  contract.summaryColumns.texture_summary));
assert(isequal(textureSummary.weighting, ...
  ["pixel_weighted"; "area_weighted_grain_mean"]));
assert(isequal(string(cAxisDistribution.Properties.VariableNames), ...
  contract.summaryColumns.c_axis_orientation_distribution));
assert(height(cAxisDistribution) == 12);
assert(all(cAxisDistribution.c_axis_ad_acute_deg >= 0 & ...
  cAxisDistribution.c_axis_ad_acute_deg <= 90));
assert(all(isfinite(cAxisDistribution.c_axis_azimuth_about_ad_deg)));
assert(max(abs(cAxisDistribution.c_axis_x_ad.^2 + ...
  cAxisDistribution.c_axis_y_td_rd.^2 + ...
  cAxisDistribution.c_axis_z_nd.^2 - 1)) < 1e-12);
canonicalRows = cAxisDistribution.weighting == "pixel_weighted";
assert(nnz(canonicalRows) == 6);
canonicalRows = cAxisDistribution(canonicalRows,:);
expectedComponents = [1 0 0; 1 0 0; 0 1 0; 0 1 0; 0 0 1; 0 0 1];
actualComponents = [canonicalRows.c_axis_x_ad, ...
  canonicalRows.c_axis_y_td_rd, canonicalRows.c_axis_z_nd];
assert(max(abs(actualComponents - expectedComponents), [], "all") < 1e-12);
assert(max(abs(canonicalRows.c_axis_ad_acute_deg - ...
  [0; 0; 90; 90; 90; 90])) < 1e-12);
assert(max(abs(canonicalRows.c_axis_azimuth_about_ad_deg - ...
  [0; 0; 0; 0; 90; 90])) < 1e-12);

% Independent c-axis mixture oracle. Pixel orientations have 1:1 AD:TD
% weights, while the same two grain means have areas 1:3.
mixtureOptions = textureOptions;
mixtureOptions.pixel_ids = [1; 2];
mixtureOptions.grain_ids = [1; 2];
[mixtureSummary, mixtureDistribution, mixtureModel] = ...
  compute_texture_metrics([cAlongAd; cAlongTd], ...
  [cAlongAd; cAlongTd], [1; 3], meta, mixtureOptions);
pixelMixture = mixtureSummary.weighting == "pixel_weighted";
grainMixture = mixtureSummary.weighting == "area_weighted_grain_mean";
assert(abs(mixtureSummary.c_axis_ad_mean_deg(pixelMixture) - 45) < 1e-12);
assert(abs(mixtureSummary.c_axis_ad_mean_deg(grainMixture) - 67.5) < 1e-12);
assert(abs(mixtureSummary.c_axis_within_15deg_fraction(grainMixture) - .25) ...
  < 1e-12);
assert(abs(mixtureSummary.c_axis_within_30deg_fraction(grainMixture) - .25) ...
  < 1e-12);
assert(abs(mixtureSummary.c_axis_ad_median_deg(grainMixture) - 90) < 1e-12);
grainMixtureRows = mixtureDistribution.weighting == ...
  "area_weighted_grain_mean";
assert(isequal(mixtureDistribution.source_weight(grainMixtureRows), [1; 3]));
assert(double(norm(mixtureModel.pixel_weighted_odf - ...
  mixtureModel.area_weighted_grain_mean_odf)) > 0);
assert(isfield(textureModel, "pixel_weighted_odf"));
assert(isfield(textureModel, "area_weighted_grain_mean_odf"));
assert(isa(textureModel.pixel_weighted_odf, "SO3FunHarmonic"));
assert(isa(textureModel.area_weighted_grain_mean_odf, ...
  "SO3FunHarmonic"));
expectedBandwidth = SO3DeLaValleePoussinKernel("halfwidth", ...
  contract.parameters.texture_kernel_halfwidth_deg * degree).bandwidth;
assert(textureModel.harmonic_bandwidth == expectedBandwidth);
assert(isfinite(textureModel.pixel_weighted_density_normalization_factor) && ...
  textureModel.pixel_weighted_density_normalization_factor > 0);
assert(isfinite( ...
  textureModel.area_weighted_grain_mean_density_normalization_factor) && ...
  textureModel.area_weighted_grain_mean_density_normalization_factor > 0);
assert(all(isfinite(textureSummary{:, 8:end}), "all"));
assert(isreal(textureSummary.clipped_grid_entropy));
assert(isfield(textureModel, "pixel_weighted_entropy_diagnostics"));
assert(isfield(textureModel, ...
  "area_weighted_grain_mean_entropy_diagnostics"));
diagnostics = textureModel.pixel_weighted_entropy_diagnostics;
assert(diagnostics.grid_point_count > 0);
assert(diagnostics.negative_grid_fraction >= 0 && ...
  diagnostics.negative_grid_fraction <= 1);
assert(diagnostics.clipped_negative_l1_fraction >= 0 && ...
  diagnostics.clipped_negative_l1_fraction <= 1);

% A small deterministic 24-orientation fixture avoids stochastic tests. It
% is only a synthetic invariant set; production uses every valid orientation.
[phi1, Phi, phi2] = ndgrid((0:5) * 60 * degree, ...
  [22.5 67.5] * degree, [0 90] * degree);
uniformOrientations = orientation.byEuler(phi1(:), Phi(:), phi2(:), cs);
singleOrientations = repmat(orientation.id(cs), 24, 1);
randomTextureOptions = rmfield(textureOptions, {'pixel_ids', 'grain_ids'});

uniformMeta = meta;
uniformMeta.sample = "uniform";
singleMeta = meta;
singleMeta.sample = "single";
[uniformSummary, ~] = compute_texture_metrics(uniformOrientations, ...
  uniformOrientations, (1:24)', uniformMeta, randomTextureOptions);
[singleSummary, ~] = compute_texture_metrics(singleOrientations, ...
  singleOrientations, ones(24, 1), singleMeta, randomTextureOptions);
assert(uniformSummary.texture_index(1) < singleSummary.texture_index(1));
assert(uniformSummary.max_mrd(1) < singleSummary.max_mrd(1));
assert(uniformSummary.clipped_grid_entropy(1) > ...
  singleSummary.clipped_grid_entropy(1));
assert(all(isfinite(uniformSummary{:, 8:end}), "all"));
assert(all(isfinite(singleSummary{:, 8:end}), "all"));
assert(isreal(uniformSummary.clipped_grid_entropy) && ...
  isreal(singleSummary.clipped_grid_entropy));

% Focused regression for MTEX's absolute weight-scale behavior: the
% unit-mean density weights must make an ODF invariant to multiplying every
% physical grain area by the same constant. Avoid repeating M-index and the
% full entropy-grid evaluation merely to check the construction invariant.
focusedKernel = SO3DeLaValleePoussinKernel("halfwidth", ...
  contract.parameters.texture_kernel_halfwidth_deg * degree);
physicalWeights = (1:24)';
scaledWeights = 10 * physicalWeights;
densityWeights = physicalWeights / mean(physicalWeights);
scaledDensityWeights = scaledWeights / mean(scaledWeights);
focusedOdf = calcDensity(uniformOrientations, "kernel", focusedKernel, ...
  "weights", densityWeights, "silent");
scaledFocusedOdf = calcDensity(uniformOrientations, "kernel", ...
  focusedKernel, "weights", scaledDensityWeights, "silent");
focusedOdf = SO3FunHarmonic(focusedOdf, "bandwidth", ...
  focusedKernel.bandwidth);
scaledFocusedOdf = SO3FunHarmonic(scaledFocusedOdf, "bandwidth", ...
  focusedKernel.bandwidth);
assert(abs(double(mean(focusedOdf)) - 1) < 1e-6);
assert(abs(double(mean(scaledFocusedOdf)) - 1) < 1e-6);
assert(double(norm(focusedOdf - scaledFocusedOdf)) < 1e-12);
assert(abs(double(norm(focusedOdf)^2) - ...
  double(norm(scaledFocusedOdf)^2)) < 1e-12);

if scanRoot ~= ""
  assert(outputDir ~= "", "outputDir is required for integration testing.");
  assert(~isempty(which("generate_comprehensive_intragranular_texture")), ...
    "Task 4 generator is missing.");
  pairOutput = fullfile(outputDir, "raw_denoised_pair");
  [pairIntra, pairTexture, pairDistance] = ...
    generate_comprehensive_intragranular_texture(scanRoot, pairOutput, ...
    struct("catalog_rows", [1 2]));
  assert(height(pairIntra) == 8);
  assert(height(pairTexture) == 4);
  assert(height(pairDistance) == 2);
  assert(all(isfinite(pairDistance.odf_harmonic_l2_distance)));
  verify_artifacts(pairOutput, contract, 2);
  verify_grain_texture_scope(pairOutput, contract);

  recoveryOptions = struct("catalog_rows", [1 2], ...
    "render_odf_only", true, ...
    "texture_render_output_name", "odf_sections.recovered.png");
  [recoveryIntra, recoveryTexture, recoveryDistance] = ...
    generate_comprehensive_intragranular_texture(scanRoot, pairOutput, ...
    recoveryOptions);
  assert(isempty(recoveryIntra) && isempty(recoveryTexture) && ...
    isempty(recoveryDistance));
  recoveredOdfPath = fullfile(pairOutput, contract.directories(6), ...
    recoveryOptions.texture_render_output_name);
  recoveredOdfInfo = dir(recoveredOdfPath);
  assert(isscalar(recoveredOdfInfo) && recoveredOdfInfo.bytes > 0);

  [intraAll, textureAll, odfDistance] = ...
    generate_comprehensive_intragranular_texture(scanRoot, outputDir);
  assert(height(intraAll) == 48);
  assert(height(textureAll) == 24);
  assert(height(odfDistance) == 12);
  assert(all(isfinite(intraAll{:, 9:end}), "all"));
  assert(all(isfinite(textureAll{:, 9:end}), "all"));
  assert(all(isfinite(odfDistance.odf_harmonic_l2_distance)));
  assert(all(odfDistance.kernel_halfwidth_deg == 5));
  verify_artifacts(outputDir, contract, 12);
end

fprintf("test_comprehensive_intragranular_texture passed\n");
end

function verify_grain_texture_scope(outputDir, contract)
grainPath = fullfile(outputDir, contract.directories(5), ...
  "intragranular_by_grain.csv");
cAxisPath = fullfile(outputDir, contract.directories(6), ...
  "c_axis_orientation_distribution.csv");
grainOptions = detectImportOptions(grainPath);
grainOptions.SelectedVariableNames = ["sample", "variant", "grain_id", ...
  "num_pixels"];
cAxisOptions = detectImportOptions(cAxisPath);
cAxisOptions.SelectedVariableNames = ["sample", "variant", "weighting", ...
  "source_id"];
grainRows = readtable(grainPath, grainOptions);
cAxisRows = readtable(cAxisPath, cAxisOptions);
grainRows.sample = string(grainRows.sample);
grainRows.variant = string(grainRows.variant);
cAxisRows.sample = string(cAxisRows.sample);
cAxisRows.variant = string(cAxisRows.variant);
cAxisRows.weighting = string(cAxisRows.weighting);
for sampleName = unique(grainRows.sample, "stable")'
  for variant = ["raw", "denoised"]
    grainCondition = grainRows.sample == sampleName & ...
      grainRows.variant == variant;
    textureCondition = cAxisRows.sample == sampleName & ...
      cAxisRows.variant == variant & ...
      cAxisRows.weighting == "area_weighted_grain_mean";
    expectedIds = sort(double(grainRows.grain_id(grainCondition)));
    actualIds = sort(double(cAxisRows.source_id(textureCondition)));
    assert(~isempty(expectedIds) && numel(unique(expectedIds)) == ...
      numel(expectedIds));
    assert(all(grainRows.num_pixels(grainCondition) >= ...
      contract.parameters.min_grain_pixels));
    assert(isequal(actualIds, expectedIds), ...
      ["Area-weighted grain texture includes grains outside the " ...
      "registered min_grain_pixels scope."]);
  end
end
end

function verify_artifacts(outputDir, contract, expectedDistanceRows)
intraDir = fullfile(outputDir, contract.directories(5));
textureDir = fullfile(outputDir, contract.directories(6));
for fileName = contract.artifacts.dir_04_intragranular'
  assert(isfile(fullfile(intraDir, fileName)), ...
    "Missing intragranular artifact: %s", fileName);
end
for fileName = contract.artifacts.dir_05_texture'
  assert(isfile(fullfile(textureDir, fileName)), ...
    "Missing texture artifact: %s", fileName);
end
distancePath = fullfile(textureDir, "raw_denoised_odf_distance.csv");
assert(isfile(distancePath));
diagnosticPath = fullfile(textureDir, ...
  "texture_numerical_diagnostics.csv");
assert(isfile(diagnosticPath));
diagnosticTable = readtable(diagnosticPath, "TextType", "string");
assert(isequal(string(diagnosticTable.Properties.VariableNames), ...
  contract.summaryColumns.texture_numerical_diagnostics));
assert(all(diagnosticTable.total_ti_hex_grain_count == ...
  diagnosticTable.retained_ti_hex_grain_count + ...
  diagnosticTable.excluded_small_grain_count));
assert(all(diagnosticTable.total_ti_hex_grain_count >= 1 & ...
  diagnosticTable.total_ti_hex_grain_count == ...
  fix(diagnosticTable.total_ti_hex_grain_count)));
assert(all(diagnosticTable.retained_ti_hex_grain_count >= 1 & ...
  diagnosticTable.retained_ti_hex_grain_count == ...
  fix(diagnosticTable.retained_ti_hex_grain_count)));
assert(all(diagnosticTable.excluded_small_grain_count >= 0 & ...
  diagnosticTable.excluded_small_grain_count == ...
  fix(diagnosticTable.excluded_small_grain_count)));
assert(all(diagnosticTable.excluded_small_grain_area_fraction >= 0 & ...
  diagnosticTable.excluded_small_grain_area_fraction <= 1));
assert(all(isfinite(diagnosticTable.density_normalization_factor) & ...
  diagnosticTable.density_normalization_factor > 0));
assert(all(diagnosticTable.negative_grid_fraction >= 0 & ...
  diagnosticTable.negative_grid_fraction <= 1));
assert(all(diagnosticTable.clipped_negative_l1_fraction >= 0 & ...
  diagnosticTable.clipped_negative_l1_fraction <= 1));
assert(isscalar(expectedDistanceRows) && expectedDistanceRows >= 0 && ...
  expectedDistanceRows == fix(expectedDistanceRows));
distanceTable = readtable(distancePath, "TextType", "string");
assert(isequal(string(distanceTable.Properties.VariableNames), ...
  contract.summaryColumns.raw_denoised_odf_distance));
assert(height(distanceTable) == expectedDistanceRows);
parameterPath = fullfile(textureDir, "task4_analysis_parameters.csv");
parameterTable = readtable(parameterPath, "TextType", "string");
assert(isequal(string(parameterTable.Properties.VariableNames), ...
  contract.summaryColumns.task4_analysis_parameters));
assert(height(parameterTable) == expectedDistanceRows);
conditionKeys = parameterTable.sample + "|" + parameterTable.variant;
assert(numel(unique(conditionKeys)) == height(parameterTable));
colorLimitPath = fullfile(textureDir, "texture_plot_color_limits.csv");
colorLimits = readtable(colorLimitPath, "TextType", "string");
assert(isequal(string(colorLimits.Properties.VariableNames), ...
  contract.summaryColumns.texture_plot_color_limits));
assert(height(colorLimits) == 6);
colorKeys = colorLimits.plot_kind + "|" + colorLimits.weighting;
assert(numel(unique(colorKeys)) == 6);
assert(all(colorLimits.scan_count == expectedDistanceRows));
assert(all(colorLimits.grid_resolution_deg == ...
  contract.parameters.texture_grid_resolution_deg));
assert(all(colorLimits.color_limit_min_mrd == 0));
assert(all(isfinite(colorLimits.color_limit_max_mrd) & ...
  colorLimits.color_limit_max_mrd > 0));
end
