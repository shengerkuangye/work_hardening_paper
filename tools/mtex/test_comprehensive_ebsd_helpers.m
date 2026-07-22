function test_comprehensive_ebsd_helpers(scanRoot)
%TEST_COMPREHENSIVE_EBSD_HELPERS Focused deterministic foundation tests.

arguments
  scanRoot (1,1) string
end

knownFile = fullfile(scanRoot, "d7", "ebsd_sample_7_map_15.ctf");
knownDigest = ...
  "3d7b0b60e9b23b0afa2a3073f54bb844453fd3a06923d92fa76e9abd469bf8c7";
assert(sha256_file(knownFile) == knownDigest, ...
  "MATLAB SHA-256 does not match the known PowerShell Get-FileHash digest.");

rawPhase = [1;1;2;2;0];
denoisedPhase = [1;2;2;0;0];
rawIndexed = logical([1;1;1;1;0]);
denoisedIndexed = logical([1;1;1;0;0]);
orientationCandidateDeg = [10;20;30;40;50];
rawQuality = struct("mad", (1:5).', "bc", (11:15).', ...
  "bs", (21:25).', "bands", (6:10).', "error", zeros(5,1));
denoisedQuality = rawQuality;
denoisedQuality.mad(1) = denoisedQuality.mad(1) + 0.5;
denoisedQuality.bc(2) = denoisedQuality.bc(2) + 1;
denoisedQuality.bs(3) = denoisedQuality.bs(3) - 2;
denoisedQuality.bands(4) = denoisedQuality.bands(4) + 1;
denoisedQuality.error(5) = denoisedQuality.error(5) + 4;
orientationToleranceDeg = ...
  comprehensive_ebsd_output_contract().parameters. ...
  orientation_change_tolerance_deg;

changes = compare_raw_denoised_ebsd_arrays(rawPhase, denoisedPhase, ...
  rawIndexed, denoisedIndexed, 1, 1, orientationCandidateDeg, ...
  rawQuality, denoisedQuality, orientationToleranceDeg);
assert(isequal(changes.phase_changed, logical([0;1;0;1;0])));
assert(isequal(changes.indexing_changed, logical([0;0;0;1;0])));
assert(isequal(changes.common_ti_hex, logical([1;0;0;0;0])));
assert(changes.orientation_change_deg(1) == 10);
assert(isnan(changes.orientation_change_deg(3)), ...
  "Cubic indexed pixels must not enter Ti-Hex orientation changes.");
assert(nnz(changes.orientation_changed) == 1);
assert(isequal(changes.mad_difference, [0.5;0;0;0;0]));
assert(isequal(changes.bc_difference, [0;1;0;0;0]));
assert(isequal(changes.bs_difference, [0;0;-2;0;0]));
assert(isequal(changes.bands_difference, [0;0;0;1;0]));
assert(isequal(changes.error_difference, [0;0;0;0;4]));

hexagonalTi = crystalSymmetry("6/mmm", [2.954 2.954 4.729], ...
  "mineral", "Ti-Hex");
identicalOrientation = orientation.byEuler(123.4567 * degree, ...
  89.9999 * degree, 359.9999 * degree, hexagonalTi);
selfDistanceDeg = double(angle(identicalOrientation, ...
  identicalOrientation) / degree);
assert(selfDistanceDeg > 2e-6 && selfDistanceDeg < 3e-6 && ...
  selfDistanceDeg < orientationToleranceDeg, ...
  "MTEX self-distance regression no longer exercises tolerance handling.");
smallOrientation = orientation.byEuler((123.4567 + 5e-5) * degree, ...
  89.9999 * degree, 359.9999 * degree, hexagonalTi);
changedOrientation = orientation.byEuler((123.4567 + 1e-2) * degree, ...
  89.9999 * degree, 359.9999 * degree, hexagonalTi);
smallChangeDeg = double(angle(identicalOrientation, ...
  smallOrientation) / degree);
realChangeDeg = double(angle(identicalOrientation, ...
  changedOrientation) / degree);
assert(smallChangeDeg > selfDistanceDeg && ...
  smallChangeDeg < orientationToleranceDeg);
assert(realChangeDeg > orientationToleranceDeg);
smallQuality = zero_quality(3);
toleranceChanges = compare_raw_denoised_ebsd_arrays(ones(3,1), ...
  ones(3,1), true(3,1), true(3,1), 1, 1, ...
  [selfDistanceDeg; smallChangeDeg; realChangeDeg], ...
  smallQuality, smallQuality, ...
  orientationToleranceDeg);
assert(isequal(toleranceChanges.orientation_changed, ...
  logical([0;0;1])), ...
  "Self-distance and sub-tolerance changes must not be flagged.");

[filteredOuter, removedOuter] = ...
  validate_filter_boundary_endpoint_ids([1 2;0 3;4 0], "outer");
assert(isequal(filteredOuter, [1 2]));
assert(isequal(removedOuter, logical([0;1;1])));
assert_throws_endpoint_error([0 0], "outer", ...
  "validate_filter_boundary_endpoint_ids:DegenerateOuterEndpoint");
invalidEndpoints = {[-1 2], [NaN 2], [Inf 2], [1.5 2]};
for invalidIndex = 1:numel(invalidEndpoints)
  assert_throws_endpoint_error(invalidEndpoints{invalidIndex}, "outer", ...
    "validate_filter_boundary_endpoint_ids:InvalidEndpointIds");
end
assert_throws_endpoint_error([0 2], "inner", ...
  "validate_filter_boundary_endpoint_ids:UnexpectedZeroEndpoint");

metricMinimum = [0 0 0 -2 -10 0; 0 0 0 -5 -4 -3];
metricMaximum = [1 1 4 3 6 2; 1 1 10 1 7 8];
colorLimits = comprehensive_ebsd_change_map_limits( ...
  metricMinimum, metricMaximum);
assert(isequal(colorLimits, [0 1;0 1;0 10;-5 5;-10 10;-8 8]));

fprintf("test_comprehensive_ebsd_helpers passed\n");
end

function quality = zero_quality(pointCount)
values = zeros(pointCount, 1);
quality = struct("mad", values, "bc", values, "bs", values, ...
  "bands", values, "error", values);
end

function assert_throws_endpoint_error(endpointIds, boundaryKind, ...
  expectedIdentifier)
didThrow = false;
try
  validate_filter_boundary_endpoint_ids(endpointIds, boundaryKind);
catch exception
  didThrow = exception.identifier == expectedIdentifier;
end
assert(didThrow, "Expected endpoint validation error %s.", ...
  expectedIdentifier);
end
