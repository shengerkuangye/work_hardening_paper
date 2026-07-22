function test_comprehensive_axial_propensity(scanRoot, outputRoot)
%TEST_COMPREHENSIVE_AXIAL_PROPENSITY Verify AD tensile propensity outputs.

arguments
  scanRoot (1,1) string = ...
    "data/ebsd_kpl_250221_7_df/scans"
  outputRoot (1,1) string = ...
    ".codex_tmp/task5-axial-propensity-test"
end

assert(~isempty(which("slipSystem")), ...
  "MTEX is not loaded in the current MATLAB session.");

cs = crystalSymmetry("6/mmm", [2.954 2.954 4.729], ...
  "mineral", "Ti-Hex");
orientations = [orientation.id(cs); ...
  orientation.byEuler(17 * degree, 43 * degree, 71 * degree, cs); ...
  orientation.byEuler(91 * degree, 28 * degree, 12 * degree, cs)];
grainIds = [101; 205; 309];
grainAreas = [2; 3; 5];

familyConstructors = { ...
  @slipSystem.basal, @slipSystem.prismaticA, ...
  @slipSystem.pyramidalA, @slipSystem.pyramidalCA, ...
  @slipSystem.twinT1, @slipSystem.twinC1};
for familyIndex = 1:numel(familyConstructors)
  baseSystem = familyConstructors{familyIndex}(cs);
  symSystems = baseSystem.symmetrise("antipodal");
  baseRotated = orientations * baseSystem;
  symRotated = orientations * symSystems;
  baseSchmid = reshape(abs(baseRotated.SchmidFactor(vector3d.X)), ...
    numel(orientations), []);
  symSchmid = reshape(abs(symRotated.SchmidFactor(vector3d.X)), ...
    numel(orientations), []);
  assert(all(symSchmid >= 0 & symSchmid <= 0.5 + 1e-12, "all"));
  assert(all(max(symSchmid, [], 2) + 1e-12 >= ...
    max(baseSchmid, [], 2)));
end

assert_validated_twin_definitions(cs);

options = struct("include_twins", true, "taylor_crss", []);
[byGrain, summary, meta] = compute_axial_propensity_metrics( ...
  orientations, grainIds, grainAreas, options);
expectedFamilies = ["basal_a"; "prismatic_a"; "pyramidal_a"; ...
  "pyramidal_ca"; "extension_twin_t1"; "contraction_twin_c1"];
expectedByGrainColumns = ["grain_id", "family", "max_abs_schmid", ...
  "winning_variant", "area_um2", "area_weight", "phi1_deg", ...
  "Phi_deg", "phi2_deg", "taylor_factor"];
expectedSummaryColumns = ["family", "assumption", "grain_count", ...
  "area_weighted_mean_max_schmid", ...
  "area_weighted_median_max_schmid", ...
  "area_fraction_schmid_ge_0_4", ...
  "taylor_factor_area_weighted_mean"];
assert(isequal(string(byGrain.Properties.VariableNames), ...
  expectedByGrainColumns));
assert(isequal(string(summary.Properties.VariableNames), ...
  expectedSummaryColumns));
assert(height(byGrain) == numel(grainIds) * numel(expectedFamilies));
assert(isequal(summary.family, expectedFamilies));
assert(meta.load_axis == "AD" && meta.load_axis_vector == "vector3d.X");
assert(meta.taylor_status == "not_requested_no_crss");
assert(all(isnan(byGrain.taylor_factor)));
assert(all(isnan(summary.taylor_factor_area_weighted_mean)));

for familyIndex = 1:numel(expectedFamilies)
  rows = byGrain.family == expectedFamilies(familyIndex);
  assert(isequal(byGrain.grain_id(rows), grainIds));
  assert(abs(sum(byGrain.area_weight(rows)) - 1) < 1e-12);
  assert(isequal(byGrain.area_weight(rows), grainAreas / sum(grainAreas)));
  assert(all(byGrain.max_abs_schmid(rows) >= 0 & ...
    byGrain.max_abs_schmid(rows) <= 0.5 + 1e-12));
  assert(all(byGrain.winning_variant(rows) >= 1 & ...
    byGrain.winning_variant(rows) == ...
    fix(byGrain.winning_variant(rows))));
end

taylorOptions = struct("include_twins", false, ...
  "taylor_crss", [1 1 1 1]);
[taylorByGrain, taylorSummary, taylorMeta] = ...
  compute_axial_propensity_metrics(orientations(1:2), ...
  grainIds(1:2), grainAreas(1:2), taylorOptions);
assert(taylorMeta.taylor_status == "ok_equal_crss");
assert(all(isfinite(taylorByGrain.taylor_factor) & ...
  taylorByGrain.taylor_factor > 0));
assert(all(isfinite(taylorSummary.taylor_factor_area_weighted_mean)));
assert(all(contains(taylorSummary.assumption, ...
  "Taylor sensitivity equal CRSS")));

assert(isfolder(scanRoot));
[fullByGrain, fullSummary] = generate_comprehensive_axial_propensity( ...
  scanRoot, outputRoot);
contract = comprehensive_ebsd_output_contract();
assert(isequal(string(fullSummary.Properties.VariableNames), ...
  contract.summaryColumns.axial_propensity_summary));
assert(height(fullSummary) == 12 * numel(expectedFamilies));
assert(all(fullSummary.grain_count > 0));
assert(isequal(unique(fullSummary.variant, "stable"), ...
  ["raw"; "denoised"]));
assert(all(isnan(fullSummary.taylor_factor_area_weighted_mean)));

scanKeys = unique(fullSummary(:, ["sample", "variant"]), "rows");
assert(height(scanKeys) == 12);
for scanIndex = 1:height(scanKeys)
  scanRows = fullSummary.sample == scanKeys.sample(scanIndex) & ...
    fullSummary.variant == scanKeys.variant(scanIndex);
  assert(isequal(fullSummary.family(scanRows), expectedFamilies));
  grainRows = fullByGrain.sample == scanKeys.sample(scanIndex) & ...
    fullByGrain.variant == scanKeys.variant(scanIndex);
  for familyIndex = 1:numel(expectedFamilies)
    familyRows = grainRows & ...
      fullByGrain.family == expectedFamilies(familyIndex);
    assert(nnz(familyRows) == fullSummary.grain_count( ...
      scanRows & fullSummary.family == expectedFamilies(familyIndex)));
    assert(abs(sum(fullByGrain.area_weight(familyRows)) - 1) < 1e-10);
  end
end

for sampleName = unique(fullSummary.sample, "stable")'
  pairRows = fullSummary.sample == sampleName;
  assert(nnz(fullSummary.variant(pairRows) == "raw") == ...
    numel(expectedFamilies));
  assert(nnz(fullSummary.variant(pairRows) == "denoised") == ...
    numel(expectedFamilies));
end

artifactDir = fullfile(outputRoot, "06_axial_propensity");
assert(isfile(fullfile(artifactDir, "axial_propensity_by_grain.csv")));
assert(isfile(fullfile(artifactDir, "axial_propensity_summary.csv")));
assert(isfile(fullfile(artifactDir, "axial_propensity_trends.png")));

fprintf("test_comprehensive_axial_propensity passed\n");
end

function assert_validated_twin_definitions(cs)
extensionTwin = slipSystem.twinT1(cs);
contractionTwin = slipSystem.twinC1(cs);
assert(isequal(round(extensionTwin.b.UVTW), [1 -1 0 1]));
assert(isequal(round(extensionTwin.n.hkil), [-1 1 0 2]));
assert(isequal(round(contractionTwin.b.UVTW), [-1 1 0 -2]));
assert(isequal(round(contractionTwin.n.hkil), [-1 1 0 1]));
assert(abs(dot(extensionTwin.b, extensionTwin.n, ...
  "noSymmetry")) < 1e-12);
assert(abs(dot(contractionTwin.b, contractionTwin.n, ...
  "noSymmetry")) < 1e-12);
end
