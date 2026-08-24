function test_generate_ipf_xyz_gb_maps(scanRoot)
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")),"MTEX must be loaded.");
outputDir = string(tempname);
mkdir(outputDir);
cleanupOutput = onCleanup(@() rmdir(outputDir,"s"));

summary = generate_ipf_xyz_gb_maps(scanRoot,outputDir);
expectedSamples = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"];
expectedReduction = [0;14.31;26.04;36;43.75;48.98];
assert(height(summary) == 6);
assert(isequal(summary.sample,expectedSamples));
assert(isequal(summary.cold_reduction_percent,expectedReduction));
assert(all(summary.variant == "raw"));
assert(all(summary.valid_ti_hex_point_count > 0));
assert(all(summary.grain_detection_deg == 2));
assert(all(summary.boundary_classification_deg == 10));
assert(all(summary.lagb_2_10_length_um > 0));
assert(all(summary.hagb_ge10_length_um > 0));
assert(all(abs(summary.lagb_2_10_length_fraction + ...
  summary.hagb_ge10_length_fraction - 1) < 1e-10));
assert(all(summary.coordinate_x == "AD"));
assert(all(summary.coordinate_y == "TD/RD"));
assert(all(summary.coordinate_z == "ND"));

for sampleIndex = 1:height(summary)
  imagePath = fullfile(outputDir,summary.sample(sampleIndex) + ...
    "_ipf_xyz_gb.png");
  assert(isfile(imagePath));
  imageInfo = imfinfo(imagePath);
  assert(imageInfo.Width > imageInfo.Height);
  [xDpi,yDpi] = image_resolution_dpi(imageInfo);
  assert(abs(xDpi - 600) < 1 && abs(yDpi - 600) < 1);
end

directionSuffixes = ["x_ad","y_td_rd","z_nd"];
for sampleIndex = 1:height(summary)
  for directionIndex = 1:numel(directionSuffixes)
    directionPath = fullfile(outputDir,summary.sample(sampleIndex) + ...
      "_ipf_" + directionSuffixes(directionIndex) + "_gb.png");
    assert(isfile(directionPath), ...
      "Missing standalone IPF map with legend: %s",directionPath);
    directionInfo = imfinfo(directionPath);
    assert(directionInfo.Width > directionInfo.Height);
    [xDpi,yDpi] = image_resolution_dpi(directionInfo);
    assert(abs(xDpi - 600) < 1 && abs(yDpi - 600) < 1);
  end
end

matrixPng = fullfile(outputDir,"ipf_xyz_gb_raw_matrix.png");
matrixPdf = fullfile(outputDir,"ipf_xyz_gb_raw_matrix.pdf");
summaryCsv = fullfile(outputDir,"ipf_xyz_gb_summary.csv");
assert(isfile(matrixPng) && isfile(matrixPdf) && isfile(summaryCsv));
matrixInfo = imfinfo(matrixPng);
assert(matrixInfo.Height > matrixInfo.Width && matrixInfo.Width >= 5000);
[xDpi,yDpi] = image_resolution_dpi(matrixInfo);
assert(abs(xDpi - 600) < 1 && abs(yDpi - 600) < 1);
assert(startsWith(fileread(matrixPdf),"%PDF-"));
assert(height(readtable(summaryCsv,"TextType","string")) == 6);

axialPng = fullfile(outputDir,"ipf_x_ad_six_state_montage.png");
axialPdf = fullfile(outputDir,"ipf_x_ad_six_state_montage.pdf");
assert(isfile(axialPng) && isfile(axialPdf));
axialInfo = imfinfo(axialPng);
assert(axialInfo.Width > axialInfo.Height && axialInfo.Width >= 6000);
[xDpi,yDpi] = image_resolution_dpi(axialInfo);
assert(abs(xDpi - 600) < 1 && abs(yDpi - 600) < 1);
assert(startsWith(fileread(axialPdf),"%PDF-"));

clear cleanupOutput
fprintf("test_generate_ipf_xyz_gb_maps passed\n");
end

function [xDpi,yDpi] = image_resolution_dpi(imageInfo)
if strcmpi(imageInfo.ResolutionUnit,"meter")
  xDpi = imageInfo.XResolution * 0.0254;
  yDpi = imageInfo.YResolution * 0.0254;
else
  assert(strcmpi(imageInfo.ResolutionUnit,"inch"));
  xDpi = imageInfo.XResolution;
  yDpi = imageInfo.YResolution;
end
end
