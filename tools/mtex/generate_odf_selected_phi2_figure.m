function [sampleSummary,peakSummary] = ...
  generate_odf_selected_phi2_figure(scanRoot,outputRoot)
%GENERATE_ODF_SELECTED_PHI2_FIGURE Render selected alpha-Ti ODF sections.

arguments
  scanRoot (1,1) string
  outputRoot (1,1) string
end
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")), "MTEX must be loaded.");
if ~isfolder(outputRoot), mkdir(outputRoot); end
[sampleSummary,peakSummary,odfs,catalog] = ...
  calculate_selected_phi2_odf_data(scanRoot);
writetable(sampleSummary, ...
  fullfile(outputRoot,"odf_selected_phi2_summary.csv"));
writetable(peakSummary, ...
  fullfile(outputRoot,"odf_selected_phi2_peak_positions.csv"));
selectedPhi2Deg = [0;30;60];
pngPath = fullfile(outputRoot,"odf_selected_phi2_sections.png");
pdfPath = fullfile(outputRoot,"odf_selected_phi2_sections.pdf");
render_selected_phi2_matrix(odfs,catalog,selectedPhi2Deg, ...
  sampleSummary.odf_maximum_mrd, ...
  sampleSummary.global_color_limit_max_mrd(1),pngPath,pdfPath);
end

function render_selected_phi2_matrix(odfs,catalog,selectedPhi2Deg, ...
  odfMaximumMrd,globalMaximumMrd,pngPath,pdfPath)
sampleCount = height(catalog);
sectionCount = numel(selectedPhi2Deg);
assert(sampleCount == 6 && sectionCount == 3);
tempRoot = string(tempname);
mkdir(tempRoot);
cleanupTemp = onCleanup(@() remove_temp_directory(tempRoot));
rowPaths = strings(sampleCount,1);
rowSizes = zeros(sampleCount,2);
renderedSectionCount = 0;

for scanIndex = 1:sampleCount
  rowPaths(scanIndex) = fullfile(tempRoot, ...
    sprintf("selected_phi2_row_%d.png",scanIndex));
  rowFigure = figure("Visible","off","Color","white", ...
    "Units","pixels","Position",[50 50 1100 520]);
  cleanupRow = onCleanup(@() close(rowFigure));
  plotSection(odfs{scanIndex},"phi2",selectedPhi2Deg * degree, ...
    "contourf","silent","layout",[1 sectionCount], ...
    "resolution",5 * degree,"colorRange",[0 globalMaximumMrd]);
  sectionAxes = findall(rowFigure,"Type","axes");
  sectionAxes = sectionAxes(arrayfun(@(axisHandle) ...
    isappdata(axisHandle,"sphericalPlot"),sectionAxes));
  assert(numel(sectionAxes) == sectionCount);
  assert(all(arrayfun(@(axisHandle) ...
    ~isempty(axisHandle.Children),sectionAxes)));
  clip_negative_contour_zdata(sectionAxes);
  for sectionIndex = 1:sectionCount
    sphericalPlotHandle = getappdata( ...
      sectionAxes(sectionIndex),"sphericalPlot");
    sphericalPlotHandle.TR.String = "";
  end
  mtexColorMap parula
  set(rowFigure,"Units","pixels","Position",[50 50 1100 520]);
  drawnow;
  exportgraphics(rowFigure,char(rowPaths(scanIndex)), ...
    "Resolution",600,"BackgroundColor","white");
  rowInfo = imfinfo(rowPaths(scanIndex));
  [rowDpiX,rowDpiY] = image_resolution_dpi(rowInfo);
  assert(abs(rowDpiX - 600) < 1 && abs(rowDpiY - 600) < 1);
  rowSizes(scanIndex,:) = [rowInfo.Width,rowInfo.Height];
  renderedSectionCount = renderedSectionCount + numel(sectionAxes);
  clear cleanupRow
end
assert(renderedSectionCount == 18);
assert(size(unique(rowSizes,"rows"),1) == 1);
crop_row_images_to_common_margins(rowPaths,36);

rowImages = cellfun(@imread,cellstr(rowPaths),"UniformOutput",false);
rowWidths = cellfun(@(imageData) size(imageData,2),rowImages);
rowHeights = cellfun(@(imageData) size(imageData,1),rowImages);
assert(numel(unique(rowWidths)) == 1 && numel(unique(rowHeights)) == 1);
targetWidth = rowWidths(1);
rowHeight = rowHeights(1);
outputDpi = 600;
labelWidthPx = max(520,round(0.20 * targetWidth));
rowGapPx = max(12,round(0.025 * rowHeight));
titleHeightPx = max(90,round(0.14 * rowHeight));
columnTitleHeightPx = max(80,round(0.12 * rowHeight));
topGapPx = max(16,round(0.025 * rowHeight));
bottomPadPx = max(24,round(0.035 * rowHeight));
colorbarGapPx = max(50,round(0.014 * targetWidth));
colorbarWidthPx = max(70,round(0.018 * targetWidth));
colorbarLabelWidthPx = max(180,round(0.055 * targetWidth));
rightPadPx = max(30,round(0.008 * targetWidth));
rowsHeightPx = sampleCount * rowHeight + ...
  (sampleCount - 1) * rowGapPx;
topPadPx = titleHeightPx + columnTitleHeightPx + topGapPx;
canvasWidthPx = labelWidthPx + targetWidth + colorbarGapPx + ...
  colorbarWidthPx + colorbarLabelWidthPx + rightPadPx;
canvasHeightPx = bottomPadPx + rowsHeightPx + topPadPx;

matrixFigure = figure("Visible","off","Color","white", ...
  "Units","inches","Position",[0.25 0.25 ...
  canvasWidthPx / outputDpi,canvasHeightPx / outputDpi]);
cleanupMatrix = onCleanup(@() close(matrixFigure));
for scanIndex = 1:sampleCount
  rowBottomPx = bottomPadPx + ...
    (sampleCount - scanIndex) * (rowHeight + rowGapPx);
  rowAxes = axes(matrixFigure,"Units","normalized","Position",[ ...
    labelWidthPx / canvasWidthPx,rowBottomPx / canvasHeightPx, ...
    targetWidth / canvasWidthPx,rowHeight / canvasHeightPx]);
  image(rowAxes,rowImages{scanIndex});
  axis(rowAxes,"off");
  set(rowAxes,"DataAspectRatio",[1 1 1], ...
    "PlotBoxAspectRatio",[targetWidth rowHeight 1], ...
    "PositionConstraint","innerposition");
  annotation(matrixFigure,"textbox",[ ...
    0,rowBottomPx / canvasHeightPx, ...
    0.96 * labelWidthPx / canvasWidthPx,rowHeight / canvasHeightPx], ...
    "String",sprintf("%.2f mm\n%.2f%% reduction\nMax ODF = %.2f MRD", ...
    catalog.diameter_mm(scanIndex), ...
    catalog.cold_reduction_percent(scanIndex), ...
    odfMaximumMrd(scanIndex)), ...
    "HorizontalAlignment","right","VerticalAlignment","middle", ...
    "Interpreter","none","FontSize",8,"Color","black", ...
    "LineStyle","none","Margin",0);
end

rowsTopPx = bottomPadPx + rowsHeightPx;
columnWidthPx = targetWidth / sectionCount;
for sectionIndex = 1:sectionCount
  columnLeftPx = labelWidthPx + ...
    (sectionIndex - 1) * columnWidthPx;
  annotation(matrixFigure,"textbox",[ ...
    columnLeftPx / canvasWidthPx, ...
    (rowsTopPx + topGapPx) / canvasHeightPx, ...
    columnWidthPx / canvasWidthPx, ...
    columnTitleHeightPx / canvasHeightPx], ...
    "String",sprintf("\\phi_2 = %d^\\circ", ...
    selectedPhi2Deg(sectionIndex)), ...
    "HorizontalAlignment","center","VerticalAlignment","middle", ...
    "Interpreter","tex","FontSize",9,"Color","black", ...
    "LineStyle","none","Margin",0);
end
annotation(matrixFigure,"textbox",[ ...
  labelWidthPx / canvasWidthPx, ...
  (rowsTopPx + topGapPx + columnTitleHeightPx) / canvasHeightPx, ...
  targetWidth / canvasWidthPx,titleHeightPx / canvasHeightPx], ...
  "String","Selected ODF sections (SS = 1)", ...
  "HorizontalAlignment","center","VerticalAlignment","middle", ...
  "Interpreter","none","FontSize",11,"Color","black", ...
  "LineStyle","none","Margin",0);

colormap(matrixFigure,parula);
colorbarLeftPx = labelWidthPx + targetWidth + colorbarGapPx;
scaleAxes = axes(matrixFigure,"Units","normalized","Position",[ ...
  colorbarLeftPx / canvasWidthPx,bottomPadPx / canvasHeightPx, ...
  1 / canvasWidthPx,rowsHeightPx / canvasHeightPx],"Visible","off");
clim(scaleAxes,[0 globalMaximumMrd]);
colorbarHandle = colorbar(scaleAxes,"eastoutside");
colorbarHandle.Units = "normalized";
colorbarHandle.Position = [ ...
  colorbarLeftPx / canvasWidthPx,bottomPadPx / canvasHeightPx, ...
  colorbarWidthPx / canvasWidthPx,rowsHeightPx / canvasHeightPx];
colorbarHandle.Label.String = "ODF intensity (MRD)";
colorbarHandle.Limits = [0 globalMaximumMrd];
colorbarHandle.FontSize = 8;
drawnow;

exportgraphics(matrixFigure,char(pngPath),"Resolution",600, ...
  "BackgroundColor","white");
exportgraphics(matrixFigure,char(pdfPath), ...
  "ContentType","image","Resolution",600, ...
  "BackgroundColor","white");
clear cleanupMatrix cleanupTemp
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

function crop_row_images_to_common_margins(paths,padding)
rowImages = cellfun(@imread,cellstr(paths),"UniformOutput",false);
referenceSize = size(rowImages{1});
assert(all(cellfun(@(imageData) ...
  isequal(size(imageData),referenceSize),rowImages)));
unionNonwhite = false(referenceSize(1),referenceSize(2));
for rowIndex = 1:numel(rowImages)
  unionNonwhite = unionNonwhite | any(rowImages{rowIndex} < 250,3);
end
[contentRows,contentColumns] = find(unionNonwhite);
assert(~isempty(contentRows));
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
