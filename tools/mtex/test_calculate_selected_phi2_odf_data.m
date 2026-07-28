function test_calculate_selected_phi2_odf_data(scanRoot)
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded for this test.");
[sampleSummary,peakSummary,odfs,catalog] = ...
  calculate_selected_phi2_odf_data(scanRoot);
assert(height(sampleSummary) == 6 && height(peakSummary) == 18);
assert(numel(odfs) == 6 && height(catalog) == 6);
assert(isequal(sampleSummary.sample, ...
  ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"]));
assert(isequal(catalog.sample,sampleSummary.sample));
assert(isequal(sampleSummary.cold_reduction_percent, ...
  [0;14.31;26.04;36;43.75;48.98]));
assert(all(sampleSummary.crystal_symmetry == "6/mmm"));
assert(all(sampleSummary.rotational_fundamental_zone == "622"));
assert(all(sampleSummary.specimen_symmetry == "1"));
assert(all(sampleSummary.odf_maximum_resolution_deg == 1));
assert(all(isfinite(sampleSummary.odf_maximum_mrd) & ...
  sampleSummary.odf_maximum_mrd > 0));
for sampleIndex = 1:6
  rows = (sampleIndex - 1) * 3 + (1:3);
  assert(isequal(peakSummary.phi2_deg(rows),[0;30;60]));
  assert(string(odfs{sampleIndex}.CS.pointGroup) == "6/mmm");
  assert(string(odfs{sampleIndex}.CS.properGroup.pointGroup) == "622");
  assert(string(odfs{sampleIndex}.SS.pointGroup) == "1");
end
assert(all(isfinite(peakSummary.section_peak_mrd) & ...
  peakSummary.section_peak_mrd > 0));
assert(all(peakSummary.section_grid_resolution_deg == 1));
globalLimit = max(sampleSummary.odf_maximum_mrd);
assert(all(sampleSummary.global_color_limit_max_mrd == globalLimit));
assert(all(peakSummary.global_color_limit_max_mrd == globalLimit));
end
