function test_generate_odf_diameter_montage(scanRoot)
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded for this test.");
originalConvention = getMTEXpref("EulerAngleConvention");
cleanupConvention = onCleanup(@() setMTEXpref( ...
  "EulerAngleConvention",originalConvention));
setMTEXpref("EulerAngleConvention","Roe");
assert(string(getMTEXpref("EulerAngleConvention")) == "Roe");
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

clipHelperPath = fullfile(fileparts(mfilename("fullpath")), ...
  "clip_negative_contour_zdata.m");
assert(isfile(clipHelperPath), ...
  "The negative contour-density clipping helper is missing.");
clipFigure = figure("Visible","off");
cleanupClipFigure = onCleanup(@() close(clipFigure));
clipAxes = axes(clipFigure);
[clipX,clipY] = meshgrid(linspace(-1,1,21));
clipZ = clipX .^ 2 + clipY .^ 2 - 0.35;
contourf(clipAxes,clipX,clipY,clipZ,[0 0.5 1 1.5]);
clipContour = findall(clipAxes,"Type","contour");
negativePointCount = nnz(clipContour.ZData < 0);
assert(negativePointCount > 0);
clippedPointCount = clip_negative_contour_zdata(clipAxes);
assert(clippedPointCount == negativePointCount);
assert(all(clipContour.ZData >= 0,"all"));
clear cleanupClipFigure

outputRoot = string(tempname);
mkdir(outputRoot);
cleanupOutput = onCleanup(@() remove_test_output(outputRoot));
[sampleSummary, sectionSummary] = generate_odf_diameter_montage(scanRoot, outputRoot);
assert(string(getMTEXpref("EulerAngleConvention")) == "Roe");
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
assert_bunge_section_maxima_match( ...
  registeredCatalog(1,:),sectionSummary(1:7,:));
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
clear cleanupConvention
assert(string(getMTEXpref("EulerAngleConvention")) == ...
  string(originalConvention));
end

function assert_bunge_section_maxima_match(catalogRow,sectionRows)
[ebsdFull,~] = load_comprehensive_ebsd_scan(catalogRow);
tiOrientations = ebsdFull("Ti-Hex").orientations;
kernel = SO3DeLaValleePoussinKernel("halfwidth",5 * degree);
rbfOdf = calcDensity(tiOrientations,"kernel",kernel, ...
  "weights",ones(numel(tiOrientations),1),"silent");
rbfOdf = normalize_positive_mean_density(rbfOdf);
odf = SO3FunHarmonic(rbfOdf,"bandwidth",kernel.bandwidth);
odf = normalize_positive_mean_density(odf);
phi1Deg = 0:5:355;
PhiDeg = 0:5:90;
[phi1Grid,PhiGrid] = meshgrid(phi1Deg,PhiDeg);
expectedSectionMaxima = zeros(7,1);
for sectionIndex = 1:7
  sectionOrientations = orientation.byEuler( ...
    phi1Grid(:) * degree,PhiGrid(:) * degree, ...
    sectionRows.phi2_deg(sectionIndex) * degree, ...
    "Bunge",tiOrientations.CS,tiOrientations.SS);
  sectionValues = real(eval(odf,sectionOrientations));
  expectedSectionMaxima(sectionIndex) = max(sectionValues);
end
assert(all(abs(sectionRows.section_maximum_mrd - ...
  expectedSectionMaxima) <= 1e-10));
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
