function [sampleSummary,sectionSummary] = ...
  generate_odf_diameter_montage(scanRoot,outputRoot)
%GENERATE_ODF_DIAMETER_MONTAGE Calculate raw alpha-Ti ODF diagnostics.

arguments
  scanRoot (1,1) string
  outputRoot (1,1) string
end

assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded.");

diameterOrder = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"];
phi2Deg = (0:10:60)';
kernelHalfwidthDeg = 5;
gridResolutionDeg = 5;

catalog = comprehensive_ebsd_catalog(scanRoot);
catalog = catalog(catalog.variant == "raw",:);
[isRegistered,catalogOrder] = ismember(diameterOrder,catalog.sample);
assert(all(isRegistered) && numel(unique(catalogOrder)) == 6);
catalog = catalog(catalogOrder,:);
assert(height(catalog) == 6);
assert(isequal(catalog.sample,diameterOrder));

if ~isfolder(outputRoot)
  mkdir(outputRoot);
end

kernel = SO3DeLaValleePoussinKernel( ...
  "halfwidth",kernelHalfwidthDeg * degree);
sampleCount = height(catalog);
sectionCount = numel(phi2Deg);
odfs = cell(sampleCount,1);
validTiHexOrientationCount = zeros(sampleCount,1);
sectionMaximumMrd = zeros(sampleCount,sectionCount);

phi1Deg = 0:gridResolutionDeg:355;
PhiDeg = 0:gridResolutionDeg:90;
[phi1Grid,PhiGrid] = meshgrid(phi1Deg,PhiDeg);

for scanIndex = 1:sampleCount
  [ebsdFull,~] = load_comprehensive_ebsd_scan(catalog(scanIndex,:));
  tiEbsd = ebsdFull("Ti-Hex");
  tiOrientations = tiEbsd.orientations;
  assert(~isempty(tiOrientations));
  assert(string(tiOrientations.CS.pointGroup) == "6/mmm");
  assert(string(tiOrientations.CS.properGroup.pointGroup) == "622");
  assert(string(tiOrientations.SS.pointGroup) == "1");
  [maxPhi1,maxPhi,maxPhi2] = fundamentalRegionEuler( ...
    tiOrientations.CS,tiOrientations.SS);
  assert(isequal([maxPhi1,maxPhi,maxPhi2],[360 90 60] * degree));

  validTiHexOrientationCount(scanIndex) = numel(tiOrientations);
  rbfOdf = calcDensity(tiOrientations,"kernel",kernel, ...
    "weights",ones(numel(tiOrientations),1),"silent");
  rbfOdf = rbfOdf / double(mean(rbfOdf));
  odf = SO3FunHarmonic(rbfOdf,"bandwidth",kernel.bandwidth);
  odf = odf / double(mean(odf));
  odfs{scanIndex} = odf;

  for sectionIndex = 1:sectionCount
    sectionOrientations = orientation.byEuler( ...
      phi1Grid(:) * degree,PhiGrid(:) * degree, ...
      phi2Deg(sectionIndex) * degree, ...
      tiOrientations.CS,tiOrientations.SS);
    sectionMrd = real(eval(odf,sectionOrientations));
    sectionMrd = sectionMrd(:);
    assert(~isempty(sectionMrd) && all(isfinite(sectionMrd)));
    sectionMaximumMrd(scanIndex,sectionIndex) = max(sectionMrd);
  end

  clear ebsdFull tiEbsd tiOrientations rbfOdf odf
end

assert(all(isfinite(sectionMaximumMrd),"all"));
assert(all(sectionMaximumMrd > 0,"all"));
globalMaximumMrd = max(sectionMaximumMrd,[],"all");
assert(isfinite(globalMaximumMrd) && globalMaximumMrd > 0);

sample = catalog.sample;
diameter_mm = catalog.diameter_mm;
cold_reduction_percent = catalog.cold_reduction_percent;
input_path = catalog.input_path;
valid_ti_hex_orientation_count = validTiHexOrientationCount;
crystal_symmetry = repmat("6/mmm",sampleCount,1);
rotational_fundamental_zone = repmat("622",sampleCount,1);
specimen_symmetry = repmat("1",sampleCount,1);
kernel_halfwidth_deg = repmat(kernelHalfwidthDeg,sampleCount,1);
grid_resolution_deg = repmat(gridResolutionDeg,sampleCount,1);
maximum_section_mrd = max(sectionMaximumMrd,[],2);
global_color_limit_max_mrd = repmat(globalMaximumMrd,sampleCount,1);
sampleSummary = table(sample,diameter_mm,cold_reduction_percent, ...
  input_path,valid_ti_hex_orientation_count,crystal_symmetry, ...
  rotational_fundamental_zone,specimen_symmetry, ...
  kernel_halfwidth_deg,grid_resolution_deg,maximum_section_mrd, ...
  global_color_limit_max_mrd);

sample = repelem(catalog.sample,sectionCount);
diameter_mm = repelem(catalog.diameter_mm,sectionCount);
cold_reduction_percent = repelem(cold_reduction_percent,sectionCount);
input_path = repelem(catalog.input_path,sectionCount);
valid_ti_hex_orientation_count = repelem( ...
  validTiHexOrientationCount,sectionCount);
crystal_symmetry = repmat("6/mmm",sampleCount * sectionCount,1);
rotational_fundamental_zone = repmat( ...
  "622",sampleCount * sectionCount,1);
specimen_symmetry = repmat("1",sampleCount * sectionCount,1);
kernel_halfwidth_deg = repmat( ...
  kernelHalfwidthDeg,sampleCount * sectionCount,1);
grid_resolution_deg = repmat( ...
  gridResolutionDeg,sampleCount * sectionCount,1);
phi2_deg = repmat(phi2Deg,sampleCount,1);
section_maximum_mrd = reshape(sectionMaximumMrd.',[],1);
global_color_limit_max_mrd = repmat( ...
  globalMaximumMrd,sampleCount * sectionCount,1);
sectionSummary = table(sample,diameter_mm,cold_reduction_percent, ...
  input_path,valid_ti_hex_orientation_count,crystal_symmetry, ...
  rotational_fundamental_zone,specimen_symmetry, ...
  kernel_halfwidth_deg,grid_resolution_deg,phi2_deg, ...
  section_maximum_mrd,global_color_limit_max_mrd);

assert(height(sampleSummary) == 6);
assert(height(sectionSummary) == 42);
assert(isequal(reshape(sectionSummary.phi2_deg,sectionCount,[]), ...
  repmat(phi2Deg,1,sampleCount)));
assert(globalMaximumMrd == max(sectionSummary.section_maximum_mrd));

writetable(sampleSummary,fullfile(outputRoot,"odf_diameter_summary.csv"));
writetable(sectionSummary, ...
  fullfile(outputRoot,"odf_diameter_section_summary.csv"));

pngPath = fullfile(outputRoot,"odf_diameter_full_sections.png");
pdfPath = fullfile(outputRoot,"odf_diameter_full_sections.pdf");
render_odf_section_matrix(odfs,catalog,phi2Deg,gridResolutionDeg, ...
  globalMaximumMrd,pngPath,pdfPath);
pngInfo = dir(pngPath);
pdfInfo = dir(pdfPath);
assert(~isempty(pngInfo) && pngInfo.bytes > 0, ...
  "ODF diagnostic PNG was not written.");
assert(~isempty(pdfInfo) && pdfInfo.bytes > 0, ...
  "ODF diagnostic PDF was not written.");
end

function render_odf_section_matrix(odfModels,catalog,phi2Deg, ...
  gridResolutionDeg,globalMaximumMrd,pngPath,pdfPath)
sampleCount = height(catalog);
sectionCount = numel(phi2Deg);
assert(sampleCount == 6 && sectionCount == 7);

tempRoot = string(tempname);
mkdir(tempRoot);
cleanupTemp = onCleanup(@() remove_temp_directory(tempRoot));
rowPaths = strings(sampleCount,1);
renderedSectionCount = 0;

for scanIndex = 1:sampleCount
  rowPaths(scanIndex) = fullfile(tempRoot, ...
    sprintf("odf_row_%d.png",scanIndex));
  rowFigure = figure("Visible","off","Color","white", ...
    "Position",[50 50 1500 520]);
  cleanupRow = onCleanup(@() close(rowFigure));
  plotSection(odfModels{scanIndex},"phi2",phi2Deg * degree, ...
    "contourf","silent","layout",[1 sectionCount], ...
    "resolution",gridResolutionDeg * degree, ...
    "colorRange",[0 globalMaximumMrd]);
  mtexColorMap parula
  set(rowFigure,"Units","pixels","Position",[50 50 1500 520]);
  drawnow;
  exportgraphics(rowFigure,char(rowPaths(scanIndex)), ...
    "Resolution",140,"BackgroundColor","white");
  clear cleanupRow
  renderedSectionCount = renderedSectionCount + sectionCount;
end
assert(renderedSectionCount == 42, ...
  "The ODF diagnostic must contain exactly 42 rendered sections.");
crop_row_images_to_common_margins(rowPaths,36);

rowImages = cellfun(@imread,cellstr(rowPaths),"UniformOutput",false);
rowWidths = cellfun(@(imageData) size(imageData,2),rowImages);
rowHeights = cellfun(@(imageData) size(imageData,1),rowImages);
targetWidth = max(rowWidths);
rowImages = cellfun(@(imageData) pad_image_to_width( ...
  imageData,targetWidth),rowImages,"UniformOutput",false);
labelPad = max(300,round(0.08 * targetWidth));

matrixFigure = figure("Visible","off","Color","white", ...
  "Units","inches","Position",[0.5 0.5 14 10]);
cleanupMatrix = onCleanup(@() close(matrixFigure));
layout = tiledlayout(matrixFigure,sampleCount,1, ...
  "Padding","compact","TileSpacing","compact");
axesHandles = gobjects(sampleCount,1);
for scanIndex = 1:sampleCount
  axesHandle = nexttile(layout);
  axesHandles(scanIndex) = axesHandle;
  rowImage = rowImages{scanIndex};
  image(axesHandle,[labelPad + 1,labelPad + targetWidth], ...
    [1,size(rowImage,1)],rowImage);
  axis(axesHandle,"image");
  axis(axesHandle,"off");
  xlim(axesHandle,[1,labelPad + targetWidth]);
  ylim(axesHandle,[1,size(rowImage,1)]);
  text(axesHandle,labelPad * 0.92,size(rowImage,1) / 2, ...
    sprintf("%.2f mm\n%.2f%% reduction", ...
    catalog.diameter_mm(scanIndex), ...
    catalog.cold_reduction_percent(scanIndex)), ...
    "HorizontalAlignment","right","VerticalAlignment","middle", ...
    "Interpreter","none","FontSize",8,"Color","black");
  clim(axesHandle,[0 globalMaximumMrd]);
end
colormap(matrixFigure,parula);
colorbarHandle = colorbar(axesHandles(end),"eastoutside");
colorbarHandle.Layout.Tile = "east";
colorbarHandle.Label.String = "ODF intensity (MRD)";
colorbarHandle.Limits = [0 globalMaximumMrd];
title(layout,sprintf( ...
  "ODF sections (SS = 1) | columns: \\phi_2 = %s deg", ...
  strjoin(string(phi2Deg.'),", ")), ...
  "Interpreter","tex","FontSize",11);
drawnow;
exportgraphics(matrixFigure,char(pngPath),"Resolution",600, ...
  "BackgroundColor","white");
exportgraphics(matrixFigure,char(pdfPath), ...
  "BackgroundColor","white");
clear cleanupMatrix cleanupTemp
end

function paddedImage = pad_image_to_width(imageData,targetWidth)
width = size(imageData,2);
assert(width <= targetWidth);
if width == targetWidth
  paddedImage = imageData;
  return
end
leftPad = floor((targetWidth - width) / 2);
rightPad = targetWidth - width - leftPad;
paddedImage = cat(2, ...
  repmat(cast(255,"like",imageData),size(imageData,1),leftPad,3), ...
  imageData, ...
  repmat(cast(255,"like",imageData),size(imageData,1),rightPad,3));
end

function crop_row_images_to_common_margins(paths,padding)
rowImages = cellfun(@imread,cellstr(paths),"UniformOutput",false);
imageSizes = cellfun(@size,rowImages,"UniformOutput",false);
referenceSize = imageSizes{1};
assert(all(cellfun(@(imageSize) isequal(imageSize,imageSizes{1}), ...
  imageSizes)),"ODF section rows must share one raster size.");
unionNonwhite = false(referenceSize(1),referenceSize(2));
for rowIndex = 1:numel(rowImages)
  imageData = rowImages{rowIndex};
  if ndims(imageData) == 2
    nonwhite = imageData < 250;
  else
    nonwhite = any(imageData < 250,3);
  end
  unionNonwhite = unionNonwhite | nonwhite;
end
[contentRows,contentColumns] = find(unionNonwhite);
assert(~isempty(contentRows),"All ODF section rows are blank.");
firstRow = max(1,min(contentRows) - padding);
lastRow = min(referenceSize(1),max(contentRows) + padding);
firstColumn = max(1,min(contentColumns) - padding);
lastColumn = min(referenceSize(2),max(contentColumns) + padding);
for rowIndex = 1:numel(rowImages)
  imageData = rowImages{rowIndex};
  imwrite(imageData(firstRow:lastRow,firstColumn:lastColumn,:), ...
    paths(rowIndex));
end
end

function remove_temp_directory(path)
if isfolder(path)
  rmdir(path,"s");
end
end
