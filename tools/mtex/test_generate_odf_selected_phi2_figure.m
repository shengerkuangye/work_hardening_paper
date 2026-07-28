function test_generate_odf_selected_phi2_figure(scanRoot)
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded for this test.");
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
[sampleSummary,peakSummary] = ...
  generate_odf_selected_phi2_figure(scanRoot,outputRoot);

sampleColumns = ["sample","diameter_mm","cold_reduction_percent", ...
  "input_path","valid_ti_hex_orientation_count","crystal_symmetry", ...
  "rotational_fundamental_zone","specimen_symmetry", ...
  "kernel_halfwidth_deg","odf_maximum_resolution_deg", ...
  "odf_maximum_mrd","phi1_max_deg","Phi_max_deg","phi2_max_deg", ...
  "global_color_limit_max_mrd"];
peakColumns = ["sample","diameter_mm","cold_reduction_percent", ...
  "phi2_deg","section_peak_mrd","phi1_peak_deg","Phi_peak_deg", ...
  "section_grid_resolution_deg","global_color_limit_max_mrd"];

assert(isequal(string(sampleSummary.Properties.VariableNames),sampleColumns));
assert(isequal(string(peakSummary.Properties.VariableNames),peakColumns));
assert(height(sampleSummary) == 6);
assert(height(peakSummary) == 18);
assert(isequal(sampleSummary.sample, ...
  ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"]));
assert(isequal(sampleSummary.cold_reduction_percent, ...
  [0;14.31;26.04;36;43.75;48.98]));
assert(all(sampleSummary.crystal_symmetry == "6/mmm"));
assert(all(sampleSummary.rotational_fundamental_zone == "622"));
assert(all(sampleSummary.specimen_symmetry == "1"));
assert(all(sampleSummary.kernel_halfwidth_deg == 5));
assert(all(sampleSummary.odf_maximum_resolution_deg == 1));
assert(all(isfinite(sampleSummary.odf_maximum_mrd) & ...
  sampleSummary.odf_maximum_mrd > 0));
assert(all(sampleSummary.phi1_max_deg >= 0 & ...
  sampleSummary.phi1_max_deg <= 360));
assert(all(sampleSummary.Phi_max_deg >= 0 & ...
  sampleSummary.Phi_max_deg <= 180));
assert(all(sampleSummary.phi2_max_deg >= 0 & ...
  sampleSummary.phi2_max_deg <= 360));

for sampleIndex = 1:6
  rows = (sampleIndex - 1) * 3 + (1:3);
  assert(isequal(peakSummary.sample(rows), ...
    repmat(sampleSummary.sample(sampleIndex),3,1)));
  assert(isequal(peakSummary.phi2_deg(rows),[0;30;60]));
end
assert(all(isfinite(peakSummary.section_peak_mrd) & ...
  peakSummary.section_peak_mrd > 0));
assert(all(peakSummary.phi1_peak_deg >= 0 & ...
  peakSummary.phi1_peak_deg < 360));
assert(all(peakSummary.Phi_peak_deg >= 0 & ...
  peakSummary.Phi_peak_deg <= 90));
assert(all(peakSummary.section_grid_resolution_deg == 1));
globalLimit = max(sampleSummary.odf_maximum_mrd);
assert(all(sampleSummary.global_color_limit_max_mrd == globalLimit));
assert(all(peakSummary.global_color_limit_max_mrd == globalLimit));

sampleCsv = fullfile(outputRoot,"odf_selected_phi2_summary.csv");
peakCsv = fullfile(outputRoot,"odf_selected_phi2_peak_positions.csv");
pngPath = fullfile(outputRoot,"odf_selected_phi2_sections.png");
pdfPath = fullfile(outputRoot,"odf_selected_phi2_sections.pdf");
assert(isfile(sampleCsv) && isfile(peakCsv));
assert(isfile(pngPath) && isfile(pdfPath));
assert(height(readtable(sampleCsv,"TextType","string")) == 6);
assert(height(readtable(peakCsv,"TextType","string")) == 18);

imageInfo = imfinfo(pngPath);
assert(imageInfo.Width > imageInfo.Height);
assert(imageInfo.Width >= 2400 && imageInfo.Height >= 1800);
[xDpi,yDpi] = image_resolution_dpi(imageInfo);
assert(abs(xDpi - 600) < 1 && abs(yDpi - 600) < 1);
imageData = imread(pngPath);
assert(mean(any(imageData < 250,3),"all") >= 0.18);

pdfInfo = dir(pdfPath);
assert(pdfInfo.bytes > 0);
pdfText = fileread(pdfPath);
assert(startsWith(pdfText,"%PDF-"));
imageTokens = regexp(pdfText, ...
  '/Subtype /Image\s+/Width ([0-9]+)\s+/Height ([0-9]+)',"tokens");
mediaBox = regexp(pdfText, ...
  '/MediaBox \[0 0 ([0-9.]+) ([0-9.]+)\]',"tokens","once");
assert(~isempty(imageTokens) && ~isempty(mediaBox));
pixelArea = sum(cellfun(@(token) ...
  str2double(token{1}) * str2double(token{2}),imageTokens));
pageArea = (str2double(mediaBox{1}) / 72) * ...
  (str2double(mediaBox{2}) / 72);
assert(sqrt(pixelArea / pageArea) >= 590);
clear cleanupOutput
end

function [xDpi,yDpi] = image_resolution_dpi(imageInfo)
if strcmpi(imageInfo.ResolutionUnit,"meter")
  xDpi = imageInfo.XResolution * 0.0254;
  yDpi = imageInfo.YResolution * 0.0254;
else
  assert(strcmpi(imageInfo.ResolutionUnit,"Inch"));
  xDpi = imageInfo.XResolution;
  yDpi = imageInfo.YResolution;
end
end

function remove_test_output(path)
if isfolder(path)
  rmdir(path,"s");
end
end
