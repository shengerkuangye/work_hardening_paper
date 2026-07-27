function test_generate_odf_diameter_montage(scanRoot)
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded for this test.");
helperPath = fullfile(fileparts(mfilename("fullpath")), ...
  "normalize_positive_mean_density.m");
assert(isfile(helperPath), ...
  "The guarded texture-density normalization helper is missing.");
[normalizedValues,normalizationFactor] = ...
  normalize_positive_mean_density([1 3]);
assert(normalizationFactor == 2);
assert(abs(mean(normalizedValues) - 1) < 1e-12);
assert_fails(@() normalize_positive_mean_density([NaN 1]), ...
  "finite and positive");
assert_fails(@() normalize_positive_mean_density([-1 -3]), ...
  "finite and positive");
assert_fails(@() normalize_positive_mean_density([1 1i]), ...
  "finite and positive");

outputRoot = string(tempname);
mkdir(outputRoot);
cleanupOutput = onCleanup(@() remove_test_output(outputRoot));
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
imageInfo = imfinfo(fullfile(outputRoot,"odf_diameter_full_sections.png"));
assert(imageInfo.Width > 0 && imageInfo.Height > 0);
assert(imageInfo.Width > imageInfo.Height, ...
  "The 6-row by 7-column SS=1 diagnostic must be horizontally oriented.");
assert(imageInfo.Width >= 2400 && imageInfo.Height >= 1800, ...
  "The diagnostic raster is smaller than the planned publication output.");
if strcmpi(imageInfo.ResolutionUnit,"meter")
  xResolutionDpi = imageInfo.XResolution * 0.0254;
  yResolutionDpi = imageInfo.YResolution * 0.0254;
else
  assert(strcmpi(imageInfo.ResolutionUnit,"Inch"));
  xResolutionDpi = imageInfo.XResolution;
  yResolutionDpi = imageInfo.YResolution;
end
assert(abs(xResolutionDpi - 600) < 1 && ...
  abs(yResolutionDpi - 600) < 1, ...
  "The diagnostic PNG must retain approximately 600 dpi metadata.");
imageData = imread(fullfile(outputRoot,"odf_diameter_full_sections.png"));
nonwhiteFraction = mean(any(imageData < 250,3),"all");
assert(nonwhiteFraction >= 0.18, ...
  "The diagnostic contains excessive white space.");
pdfInfo = dir(fullfile(outputRoot,"odf_diameter_full_sections.pdf"));
assert(pdfInfo.bytes > 0, "The diagnostic PDF is empty.");
pdfId = fopen(fullfile(outputRoot,"odf_diameter_full_sections.pdf"),"r");
assert(pdfId >= 0, "The diagnostic PDF cannot be opened.");
pdfHeader = fread(pdfId,5,"*char").';
fclose(pdfId);
assert(strcmp(pdfHeader,"%PDF-"), ...
  "The diagnostic PDF does not contain a valid PDF header.");
pdfBytes = fileread(fullfile(outputRoot,"odf_diameter_full_sections.pdf"));
imageTokens = regexp(pdfBytes, ...
  '/Subtype /Image\s+/Width ([0-9]+)\s+/Height ([0-9]+)',"tokens");
mediaBoxToken = regexp(pdfBytes, ...
  '/MediaBox \[0 0 ([0-9.]+) ([0-9.]+)\]',"tokens","once");
assert(~isempty(imageTokens) && ~isempty(mediaBoxToken), ...
  "The diagnostic PDF does not expose its raster and page dimensions.");
imagePixelArea = sum(cellfun(@(token) ...
  str2double(token{1}) * str2double(token{2}),imageTokens));
pageAreaSquareInches = (str2double(mediaBoxToken{1}) / 72) * ...
  (str2double(mediaBoxToken{2}) / 72);
pageEquivalentImageDpi = sqrt(imagePixelArea / pageAreaSquareInches);
assert(pageEquivalentImageDpi >= 590, ...
  "The diagnostic PDF contains less than 590 dpi of page raster data.");
roundTripSample = readtable(fullfile(outputRoot,"odf_diameter_summary.csv"),"TextType","string");
roundTripSection = readtable(fullfile(outputRoot,"odf_diameter_section_summary.csv"),"TextType","string");
assert(isequal(string(roundTripSample.Properties.VariableNames),sampleColumns));
assert(isequal(string(roundTripSection.Properties.VariableNames),sectionColumns));
assert(height(roundTripSample) == 6 && height(roundTripSection) == 42);
clear cleanupOutput
end

function assert_fails(operation, expectedMessage)
didFail = false;
try
  operation();
catch exception
  didFail = contains(exception.message,expectedMessage);
end
assert(didFail, ...
  "Expected operation to fail with a message containing '%s'.", ...
  expectedMessage);
end

function remove_test_output(path)
if isfolder(path)
  rmdir(path,"s");
end
end
