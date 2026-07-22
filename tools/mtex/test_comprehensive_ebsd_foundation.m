function test_comprehensive_ebsd_foundation(scanRoot)
%TEST_COMPREHENSIVE_EBSD_FOUNDATION Verify shared EBSD inputs and audit.

arguments
  scanRoot (1,1) string
end

assert(isfolder(scanRoot), "EBSD scan folder not found: %s", scanRoot);
assert(~isempty(which("EBSD")), ...
  "MTEX is not loaded in the current MATLAB session.");
contract = comprehensive_ebsd_output_contract();
catalog = comprehensive_ebsd_catalog(scanRoot);

expectedCatalogColumns = [ ...
  "sample", "diameter_mm", "cold_reduction_percent", "variant", ...
  "folder", "input_file", "input_path", "coordinate_x", ...
  "coordinate_y", "coordinate_z", "x_cells", "y_cells", ...
  "x_step_um", "y_step_um"];
assert(istable(catalog));
assert(isequal(string(catalog.Properties.VariableNames), ...
  expectedCatalogColumns));
assert(height(catalog) == 12);
assert(isequal(catalog.sample, repelem( ...
  ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"], 2)));
assert(isequal(catalog.diameter_mm, repelem( ...
  [7;6.48;6.02;5.6;5.25;5], 2)));
assert(isequal(catalog.cold_reduction_percent, repelem( ...
  [0;14.31;26.04;36.00;43.75;48.98], 2)));
assert(isequal(catalog.variant, repmat(["raw";"denoised"], 6, 1)));
assert(numel(unique(catalog.input_path)) == height(catalog));
assert(all(arrayfun(@isfile, catalog.input_path)));
assert(all(catalog.coordinate_x == contract.parameters.coordinate_x));
assert(all(catalog.coordinate_y == contract.parameters.coordinate_y));
assert(all(catalog.coordinate_z == contract.parameters.coordinate_z));
assert(all(catalog.x_cells == 600 & catalog.y_cells == 600));
assert(all(catalog.x_step_um == 0.5 & catalog.y_step_um == 0.5));

test_comprehensive_ebsd_helpers(scanRoot);

for sampleIndex = 1:6
  rawIndex = 2 * sampleIndex - 1;
  denoisedIndex = rawIndex + 1;
  [rawEbsd, rawMeta] = load_comprehensive_ebsd_scan( ...
    catalog(rawIndex, :));
  [denoisedEbsd, denoisedMeta] = load_comprehensive_ebsd_scan( ...
    catalog(denoisedIndex, :));
  assert_loader_result(rawEbsd, rawMeta);
  assert_loader_result(denoisedEbsd, denoisedMeta);
  rawXY = [double(rawEbsd.x(:)), double(rawEbsd.y(:))];
  denoisedXY = [double(denoisedEbsd.x(:)), ...
    double(denoisedEbsd.y(:))];
  assert(isequal(rawXY, denoisedXY));
  sampleMeta = catalog(rawIndex, :);
  [pairRow, pixelTable] = compare_raw_denoised_ebsd( ...
    rawEbsd, denoisedEbsd, sampleMeta);
  assert(pairRow.total_pixels == 360000);
  assert(pairRow.coordinate_mismatch_count == 0);
  assert(height(pixelTable) == pairRow.total_pixels);
  assert(all(ismember(contract.summaryColumns.raw_denoised_pair_audit, ...
    string(pairRow.Properties.VariableNames))));
  requiredPixelColumns = ["x_um", "y_um", "phase_changed", ...
    "indexing_changed", "orientation_change_deg", ...
    "orientation_changed", "mad_difference", "bc_difference", ...
    "bs_difference", "bands_difference", "error_difference"];
  assert(all(ismember(requiredPixelColumns, ...
    string(pixelTable.Properties.VariableNames))));
  assert(pairRow.phase_changed_pixels == nnz(pixelTable.phase_changed));
  assert(pairRow.indexing_changed_pixels == ...
    nnz(pixelTable.indexing_changed));
  assert(pairRow.orientation_changed_pixels == ...
    nnz(pixelTable.orientation_changed));
  assert(pairRow.mad_changed_pixels == ...
    nnz(pixelTable.mad_difference ~= 0));
  assert(pairRow.bc_changed_pixels == ...
    nnz(pixelTable.bc_difference ~= 0));
  assert(pairRow.bs_changed_pixels == ...
    nnz(pixelTable.bs_difference ~= 0));
  assert(pairRow.bands_changed_pixels == ...
    nnz(pixelTable.bands_difference ~= 0));
  assert(pairRow.error_changed_pixels == ...
    nnz(pixelTable.error_difference ~= 0));
  if sampleIndex == 1
    assert_7d_orientation_benchmark(catalog(rawIndex, :), ...
      catalog(denoisedIndex, :), pairRow, pixelTable, ...
      contract.parameters.orientation_change_tolerance_deg);
  end
  clear rawEbsd denoisedEbsd pixelTable
end

reconstructionOptions = struct("detection_threshold_deg", 2, ...
  "min_grain_pixels", 5);
[rawEbsd, ~] = load_comprehensive_ebsd_scan(catalog(1, :));
sourceIds = double(rawEbsd.id(:));
[grains, reconstructed, recon] = reconstruct_comprehensive_grains( ...
  rawEbsd, reconstructionOptions);
assert(isequal(double(reconstructed.id(:)), sourceIds));
assert(length(reconstructed) == 360000);
assert(recon.detection_threshold_deg == 2);
assert(recon.min_grain_pixels == 5);
assert(recon.all_grain_count == length(grains));
assert(recon.summary_grain_count == nnz(grains.numPixel >= 5));
assert(recon.native_grid_audit.nonlocal_endpoint_pair_count == 0);

auditDir = string(tempname);
cleanupAudit = onCleanup(@() remove_test_output(auditDir));
[inventory, pairs] = generate_comprehensive_ebsd_audit(scanRoot, auditDir);
assert(height(inventory) == 12);
assert(height(pairs) == 6);
assert(isequal(string(inventory.Properties.VariableNames), ...
  contract.summaryColumns.scan_inventory));
assert(isequal(string(pairs.Properties.VariableNames), ...
  contract.summaryColumns.raw_denoised_pair_audit));
assert(inventory.sha256(1) == ...
  "3d7b0b60e9b23b0afa2a3073f54bb844453fd3a06923d92fa76e9abd469bf8c7");
auditOutputDir = fullfile(auditDir, "00_audit");
assert(isfile(fullfile(auditOutputDir, "scan_inventory.csv")));
assert(isfile(fullfile(auditOutputDir, "raw_denoised_pair_audit.csv")));
assert(isfile(fullfile(auditOutputDir, "raw_denoised_change_maps.png")));

fprintf("test_comprehensive_ebsd_foundation passed\n");
end

function assert_7d_orientation_benchmark(rawCatalogRow, ...
  denoisedCatalogRow, pairRow, pixelTable, toleranceDeg)
rawCtf = read_ctf_numeric_rows(rawCatalogRow.input_path);
denoisedCtf = read_ctf_numeric_rows(denoisedCatalogRow.input_path);
assert(isequal(rawCtf(:, 2:3), denoisedCtf(:, 2:3)));
assert(isequal(rawCtf(:, 2:3), ...
  [pixelTable.x_um, pixelTable.y_um]), ...
  "MTEX point order differs from the native CTF row order.");

commonHex = rawCtf(:, 1) == 1 & denoisedCtf(:, 1) == 1;
rawEulerDeg = rawCtf(:, 6:8);
denoisedEulerDeg = denoisedCtf(:, 6:8);
exactSameEuler = commonHex & all(rawEulerDeg == denoisedEulerDeg, 2);
exactDifferentEuler = commonHex & any( ...
  rawEulerDeg ~= denoisedEulerDeg, 2);
assert(~any(pixelTable.orientation_changed(exactSameEuler)), ...
  "Identical CTF Euler rows must not be flagged as orientation changes.");

hexagonalTi = crystalSymmetry("6/mmm", [2.954 2.954 4.729], ...
  "mineral", "Ti-Hex");
rawOrientation = orientation.byEuler( ...
  rawEulerDeg(commonHex, 1) * degree, ...
  rawEulerDeg(commonHex, 2) * degree, ...
  rawEulerDeg(commonHex, 3) * degree, hexagonalTi);
denoisedOrientation = orientation.byEuler( ...
  denoisedEulerDeg(commonHex, 1) * degree, ...
  denoisedEulerDeg(commonHex, 2) * degree, ...
  denoisedEulerDeg(commonHex, 3) * degree, hexagonalTi);
externalChangeDeg = double(angle(rawOrientation, ...
  denoisedOrientation) / degree);
externalChangedCount = nnz(externalChangeDeg > toleranceDeg);
assert(pairRow.orientation_changed_pixels == externalChangedCount, ...
  "7d orientation summary differs from the independent CTF benchmark.");
fprintf("7d CTF orientation benchmark: exact Euler differences=%d, " + ...
  "physical changes above %.4g deg=%d\n", ...
  nnz(exactDifferentEuler), toleranceDeg, externalChangedCount);
end

function numericRows = read_ctf_numeric_rows(filePath)
numericRows = readmatrix(filePath, "FileType", "text", ...
  "Delimiter", sprintf("\t"), "NumHeaderLines", 16);
assert(size(numericRows, 1) == 360000 && size(numericRows, 2) == 11, ...
  "Unexpected CTF numeric table shape for %s.", filePath);
end

function assert_loader_result(ebsdFull, meta)
assert(length(ebsdFull) == 360000);
assert(numel(unique(double(ebsdFull.id(:)))) == 360000);
assert(meta.x_cells == 600);
assert(meta.y_cells == 600);
assert(meta.x_step_um == 0.5);
assert(meta.y_step_um == 0.5);
assert(meta.total_pixels == 360000);
assert(meta.coordinate_x == "AD");
assert(meta.coordinate_y == "TD_RD");
assert(meta.coordinate_z == "ND");
assert(isfile(meta.input_path));
assert(strlength(meta.mtex_version) > 0);
assert(meta.has_mad && meta.has_bc && meta.has_bs && ...
  meta.has_bands && meta.has_error);
end

function remove_test_output(outputDir)
if isfolder(outputDir)
  rmdir(outputDir, "s");
end
end
