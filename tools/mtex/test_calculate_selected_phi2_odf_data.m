function test_calculate_selected_phi2_odf_data(scanRoot)
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded for this test.");
originalConvention = getMTEXpref("EulerAngleConvention");
originalWorkingDirectory = pwd;
temporaryWorkingDirectory = string(tempname);
mkdir(temporaryWorkingDirectory);
cleanupState = onCleanup(@() restore_test_state( ...
  originalConvention,originalWorkingDirectory,temporaryWorkingDirectory));
setMTEXpref("EulerAngleConvention","Roe");
scanManifestBefore = recursive_file_manifest(scanRoot);
cd(temporaryWorkingDirectory);
[sampleSummary,peakSummary,odfs,catalog] = ...
  calculate_selected_phi2_odf_data(scanRoot);
temporaryContents = dir(temporaryWorkingDirectory);
temporaryContents = temporaryContents(~ismember( ...
  string({temporaryContents.name}),[".",".."]));
assert(isempty(temporaryContents), ...
  "Calculation interface must not create files or directories.");
scanManifestAfter = recursive_file_manifest(scanRoot);
assert(isequal(scanManifestAfter,scanManifestBefore), ...
  "Calculation interface must not alter the scan-root file manifest.");
cd(originalWorkingDirectory);
assert(height(sampleSummary) == 6 && height(peakSummary) == 18);
assert(numel(odfs) == 6 && height(catalog) == 6);
assert(isequal(sampleSummary.sample, ...
  ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"]));
assert(isequal(catalog.sample,sampleSummary.sample));
assert(all(catalog.variant == "raw"));
assert(numel(unique(catalog.input_path)) == 6);
assert(all(arrayfun(@isfile,catalog.input_path)));
assert(isequal(sampleSummary.cold_reduction_percent, ...
  [0;14.31;26.04;36;43.75;48.98]));
assert(all(sampleSummary.crystal_symmetry == "6/mmm"));
assert(all(sampleSummary.rotational_fundamental_zone == "622"));
assert(all(sampleSummary.specimen_symmetry == "1"));
assert(all(sampleSummary.odf_maximum_resolution_deg == 1));
assert(all(isfinite(sampleSummary.odf_maximum_mrd) & ...
  sampleSummary.odf_maximum_mrd > 0));
assert(all(sampleSummary.valid_ti_hex_orientation_count > 0));
phi1Deg = 0:359;
PhiDeg = 0:90;
[phi1Grid,PhiGrid] = meshgrid(phi1Deg,PhiDeg);
for sampleIndex = 1:6
  rows = (sampleIndex - 1) * 3 + (1:3);
  assert(isequal(peakSummary.phi2_deg(rows),[0;30;60]));
  assert(string(odfs{sampleIndex}.CS.pointGroup) == "6/mmm");
  assert(string(odfs{sampleIndex}.CS.properGroup.pointGroup) == "622");
  assert(string(odfs{sampleIndex}.SS.pointGroup) == "1");

  [expectedMaximum,expectedOrientation] = max(odfs{sampleIndex}, ...
    "resolution",1 * degree);
  assert(isscalar(expectedMaximum) && isscalar(expectedOrientation));
  assert(isfinite(expectedMaximum) && expectedMaximum > 0);
  [expectedPhi1,expectedPhi,expectedPhi2] = ...
    Euler(expectedOrientation,"Bunge");
  expectedEulerDeg = [mod(expectedPhi1 / degree,360), ...
    expectedPhi / degree,mod(expectedPhi2 / degree,360)];
  assert(all(isfinite(expectedEulerDeg)));
  assert(expectedEulerDeg(1) >= 0 && expectedEulerDeg(1) < 360);
  assert(expectedEulerDeg(2) >= 0 && expectedEulerDeg(2) <= 180);
  assert(expectedEulerDeg(3) >= 0 && expectedEulerDeg(3) < 360);
  assert(abs(sampleSummary.odf_maximum_mrd(sampleIndex) - ...
    double(expectedMaximum)) <= 1e-10);
  actualEulerDeg = [sampleSummary.phi1_max_deg(sampleIndex), ...
    sampleSummary.Phi_max_deg(sampleIndex), ...
    sampleSummary.phi2_max_deg(sampleIndex)];
  assert(all(abs(actualEulerDeg - expectedEulerDeg) <= 1e-10));

  for sectionIndex = 1:3
    sectionOrientations = orientation.byEuler( ...
      phi1Grid(:) * degree,PhiGrid(:) * degree, ...
      peakSummary.phi2_deg(rows(sectionIndex)) * degree, ...
      "Bunge",odfs{sampleIndex}.CS,odfs{sampleIndex}.SS);
    sectionValues = real(eval(odfs{sampleIndex},sectionOrientations));
    assert(iscolumn(sectionValues) && numel(sectionValues) == numel(phi1Grid));
    assert(all(isfinite(sectionValues)));
    [expectedSectionPeak,expectedPeakIndex] = max(sectionValues);
    assert(isscalar(expectedSectionPeak) && isscalar(expectedPeakIndex));
    assert(isfinite(expectedSectionPeak) && expectedSectionPeak > 0);
    expectedPhi1Peak = phi1Grid(expectedPeakIndex);
    expectedPhiPeak = PhiGrid(expectedPeakIndex);
    assert(expectedPhi1Peak >= 0 && expectedPhi1Peak <= 359);
    assert(expectedPhiPeak >= 0 && expectedPhiPeak <= 90);
    actualRow = rows(sectionIndex);
    assert(abs(peakSummary.section_peak_mrd(actualRow) - ...
      expectedSectionPeak) <= 1e-10);
    assert(peakSummary.phi1_peak_deg(actualRow) == expectedPhi1Peak);
    assert(peakSummary.Phi_peak_deg(actualRow) == expectedPhiPeak);
  end
end
assert(all(isfinite(peakSummary.section_peak_mrd) & ...
  peakSummary.section_peak_mrd > 0));
assert(all(peakSummary.section_grid_resolution_deg == 1));
globalLimit = max(sampleSummary.odf_maximum_mrd);
assert(all(sampleSummary.global_color_limit_max_mrd == globalLimit));
assert(all(peakSummary.global_color_limit_max_mrd == globalLimit));
assert(getMTEXpref("EulerAngleConvention") == "Roe");
clear cleanupState
end

function manifest = recursive_file_manifest(scanRoot)
entries = dir(fullfile(scanRoot,"**","*"));
entries = entries(~[entries.isdir]);
fullPaths = string(fullfile({entries.folder},{entries.name})).';
relative_path = extractAfter(fullPaths,strlength(scanRoot) + 1);
bytes = [entries.bytes].';
modified_datenum = [entries.datenum].';
sha256 = strings(numel(entries),1);
for fileIndex = 1:numel(entries)
  sha256(fileIndex) = sha256_file(fullPaths(fileIndex));
end
manifest = sortrows(table(relative_path,bytes,modified_datenum,sha256), ...
  "relative_path");
end

function restore_test_state(originalConvention,originalWorkingDirectory, ...
    temporaryWorkingDirectory)
setMTEXpref("EulerAngleConvention",originalConvention);
cd(originalWorkingDirectory);
if isfolder(temporaryWorkingDirectory)
  rmdir(temporaryWorkingDirectory,"s");
end
end
