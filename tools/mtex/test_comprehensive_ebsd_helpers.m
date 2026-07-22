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

changes = compare_raw_denoised_ebsd_arrays(rawPhase, denoisedPhase, ...
  rawIndexed, denoisedIndexed, 1, 1, orientationCandidateDeg, ...
  rawQuality, denoisedQuality);
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

metricMinimum = [0 0 0 -2 -10 0; 0 0 0 -5 -4 -3];
metricMaximum = [1 1 4 3 6 2; 1 1 10 1 7 8];
colorLimits = comprehensive_ebsd_change_map_limits( ...
  metricMinimum, metricMaximum);
assert(isequal(colorLimits, [0 1;0 1;0 10;-5 5;-10 10;-8 8]));

fprintf("test_comprehensive_ebsd_helpers passed\n");
end
