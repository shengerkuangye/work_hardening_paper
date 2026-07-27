function test_generate_odf_diameter_montage(scanRoot, outputRoot)
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded for this test.");
if isfolder(outputRoot), rmdir(outputRoot,"s"); end
[sampleSummary, sectionSummary] = generate_odf_diameter_montage(scanRoot, outputRoot);
sampleColumns = ["sample","diameter_mm","cold_reduction_percent", ...
  "input_path","valid_ti_hex_orientation_count","crystal_symmetry", ...
  "rotational_fundamental_zone","specimen_symmetry", ...
  "kernel_halfwidth_deg","grid_resolution_deg","maximum_section_mrd", ...
  "global_color_limit_max_mrd"];
sectionColumns = ["sample","diameter_mm","cold_reduction_percent", ...
  "input_path","valid_ti_hex_orientation_count","crystal_symmetry", ...
  "rotational_fundamental_zone","specimen_symmetry", ...
  "kernel_halfwidth_deg","grid_resolution_deg","phi2_deg", ...
  "section_maximum_mrd","global_color_limit_max_mrd"];
assert(isequal(string(sampleSummary.Properties.VariableNames),sampleColumns));
assert(isequal(string(sectionSummary.Properties.VariableNames),sectionColumns));
assert(height(sampleSummary) == 6);
assert(height(sectionSummary) == 42);
assert(isequal(sampleSummary.sample,["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"]));
assert(isequal(sampleSummary.diameter_mm,[7;6.48;6.02;5.6;5.25;5]));
rawCatalog = comprehensive_ebsd_catalog(scanRoot);
rawCatalog = rawCatalog(rawCatalog.variant == "raw",:);
[isRegistered,catalogOrder] = ismember(sampleSummary.sample,rawCatalog.sample);
assert(all(isRegistered) && numel(unique(catalogOrder)) == 6);
registeredCatalog = rawCatalog(catalogOrder,:);
assert(isequal(registeredCatalog.cold_reduction_percent, ...
  [0;14.31;26.04;36;43.75;48.98]));
assert(isequal(sampleSummary.cold_reduction_percent, ...
  registeredCatalog.cold_reduction_percent));
assert(all(sampleSummary.valid_ti_hex_orientation_count > 0));
assert(all(strlength(sampleSummary.input_path) > 0));
assert(all(arrayfun(@isfile,sampleSummary.input_path)));
for sampleIndex = 1:6
  sectionRows = (sampleIndex - 1) * 7 + (1:7);
  assert(isequal(sectionSummary.phi2_deg(sectionRows),(0:10:60)'));
  assert(isequal(sectionSummary.sample(sectionRows),repmat(sampleSummary.sample(sampleIndex),7,1)));
  assert(isequal(sectionSummary.diameter_mm(sectionRows),repmat(sampleSummary.diameter_mm(sampleIndex),7,1)));
  assert(isequal(sectionSummary.cold_reduction_percent(sectionRows),repmat(sampleSummary.cold_reduction_percent(sampleIndex),7,1)));
  assert(isequal(sectionSummary.input_path(sectionRows),repmat(sampleSummary.input_path(sampleIndex),7,1)));
  assert(isequal(sectionSummary.valid_ti_hex_orientation_count(sectionRows),repmat(sampleSummary.valid_ti_hex_orientation_count(sampleIndex),7,1)));
  assert(sampleSummary.maximum_section_mrd(sampleIndex) == max(sectionSummary.section_maximum_mrd(sectionRows)));
end
assert(all(sampleSummary.crystal_symmetry == "6/mmm"));
assert(all(sampleSummary.rotational_fundamental_zone == "622"));
assert(all(sampleSummary.specimen_symmetry == "1"));
assert(all(sampleSummary.kernel_halfwidth_deg == 5));
assert(all(sampleSummary.grid_resolution_deg == 5));
assert(all(sectionSummary.crystal_symmetry == "6/mmm"));
assert(all(sectionSummary.rotational_fundamental_zone == "622"));
assert(all(sectionSummary.specimen_symmetry == "1"));
assert(all(sectionSummary.kernel_halfwidth_deg == 5));
assert(all(sectionSummary.grid_resolution_deg == 5));
assert(all(isfinite(sectionSummary.section_maximum_mrd) & sectionSummary.section_maximum_mrd > 0));
assert(all(isfinite(sectionSummary.global_color_limit_max_mrd) & sectionSummary.global_color_limit_max_mrd > 0));
assert(numel(unique(sectionSummary.global_color_limit_max_mrd)) == 1);
assert(sectionSummary.global_color_limit_max_mrd(1) == max(sectionSummary.section_maximum_mrd));
assert(all(isfinite(sampleSummary.global_color_limit_max_mrd) & sampleSummary.global_color_limit_max_mrd > 0));
assert(all(sampleSummary.global_color_limit_max_mrd == sectionSummary.global_color_limit_max_mrd(1)));
assert(isfile(fullfile(outputRoot,"odf_diameter_summary.csv")));
assert(isfile(fullfile(outputRoot,"odf_diameter_section_summary.csv")));
assert(isfile(fullfile(outputRoot,"odf_diameter_full_sections.png")));
assert(isfile(fullfile(outputRoot,"odf_diameter_full_sections.pdf")));
roundTripSample = readtable(fullfile(outputRoot,"odf_diameter_summary.csv"),"TextType","string");
roundTripSection = readtable(fullfile(outputRoot,"odf_diameter_section_summary.csv"),"TextType","string");
assert(isequal(string(roundTripSample.Properties.VariableNames),sampleColumns));
assert(isequal(string(roundTripSection.Properties.VariableNames),sectionColumns));
assert(height(roundTripSample) == 6 && height(roundTripSection) == 42);
end
