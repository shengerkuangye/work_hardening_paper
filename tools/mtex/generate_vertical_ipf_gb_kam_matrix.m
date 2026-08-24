function summary = generate_vertical_ipf_gb_kam_matrix(scanRoot, outputDir)
%GENERATE_VERTICAL_IPF_GB_KAM_MATRIX Six-state raw-EBSD paper matrix.
% Rows follow increasing cold reduction. Columns contain IPF-X/AD,
% classified Ti-Hex boundaries, and KAM (order 1, 5 degree cutoff).

arguments
  scanRoot (1,1) string
  outputDir (1,1) string
end

assert(isfolder(scanRoot), "EBSD scan folder not found: %s", scanRoot);
assert(~isempty(which("EBSD")), "MTEX must be loaded before this function.");
if ~isfolder(outputDir)
  mkdir(outputDir);
end

catalog = comprehensive_ebsd_catalog(scanRoot);
catalog = catalog(catalog.variant == "raw", :);
[~, order] = sort(catalog.cold_reduction_percent);
catalog = catalog(order, :);
assert(height(catalog) == 6, "Expected six raw EBSD states.");

detectionDeg = 2;
hagbDeg = 15;
minGrainPixels = 5;
kamOrder = 1;
kamCutoffDeg = 5;
kamDisplayDeg = [0 5];
sampleCount = height(catalog);

xCoordinates = cell(sampleCount, 1);
yCoordinates = cell(sampleCount, 1);
ipfImages = cell(sampleCount, 1);
bcImages = cell(sampleCount, 1);
kamImages = cell(sampleCount, 1);
boundary2To5 = cell(sampleCount, 2);
boundary5To15 = cell(sampleCount, 2);
boundary15Plus = cell(sampleCount, 2);
boundaryLength2To5 = zeros(sampleCount, 1);
boundaryLength5To15 = zeros(sampleCount, 1);
boundaryLength15Plus = zeros(sampleCount, 1);
kamMean = zeros(sampleCount, 1);
kamMedian = zeros(sampleCount, 1);
kamP90 = zeros(sampleCount, 1);
bcForLimits = cell(sampleCount, 1);
tiCrystalSymmetry = [];

for sampleIndex = 1:sampleCount
  fprintf("VERTICAL_MATRIX sample=%s reduction=%.2f%%\n", ...
    catalog.sample(sampleIndex), catalog.cold_reduction_percent(sampleIndex));
  [ebsdFull, ~] = load_comprehensive_ebsd_scan(catalog(sampleIndex, :));
  tiEbsd = ebsdFull("Ti-Hex");
  assert(~isempty(tiEbsd), "Ti-Hex phase is absent.");
  if isempty(tiCrystalSymmetry)
    tiCrystalSymmetry = tiEbsd.CS;
  end

  [grains, grainId] = calcGrains(ebsdFull, "unitCell", ...
    "threshold", [hagbDeg detectionDeg] * degree, ...
    "minPixel", minGrainPixels);
  ebsdFull.grainId = grainId;
  retainedGrainIds = double(grains.id(grains.numPixel >= minGrainPixels));
  assert(~isempty(retainedGrainIds), "No reconstructed grains pass minPixel.");

  innerBoundary = grains.innerBoundary("Ti-Hex", "Ti-Hex");
  outerBoundary = grains.boundary("Ti-Hex", "Ti-Hex");
  innerAngles = double(angle(innerBoundary.misorientation) / degree);
  outerAngles = double(angle(outerBoundary.misorientation) / degree);
  boundaryA = innerBoundary(innerAngles >= 2 & innerAngles < 5);
  boundaryB = innerBoundary(innerAngles >= 5 & innerAngles < hagbDeg);
  boundaryC = outerBoundary(outerAngles >= hagbDeg);

  [boundary2To5{sampleIndex, 1}, boundary2To5{sampleIndex, 2}] = ...
    boundary_coordinates(boundaryA);
  [boundary5To15{sampleIndex, 1}, boundary5To15{sampleIndex, 2}] = ...
    boundary_coordinates(boundaryB);
  [boundary15Plus{sampleIndex, 1}, boundary15Plus{sampleIndex, 2}] = ...
    boundary_coordinates(boundaryC);
  boundaryLength2To5(sampleIndex) = boundary_length(boundaryA);
  boundaryLength5To15(sampleIndex) = boundary_length(boundaryB);
  boundaryLength15Plus(sampleIndex) = boundary_length(boundaryC);

  [xCoordinates{sampleIndex}, yCoordinates{sampleIndex}, nativeIndices] = ...
    native_grid_indices(ebsdFull);
  ipfImages{sampleIndex} = build_ipf_image(ebsdFull, tiEbsd, ...
    nativeIndices, numel(yCoordinates{sampleIndex}), ...
    numel(xCoordinates{sampleIndex}));
  bcValues = double(ebsdFull.prop.bc(:));
  bcImages{sampleIndex} = scalar_image(bcValues, nativeIndices, ...
    numel(yCoordinates{sampleIndex}), numel(xCoordinates{sampleIndex}));
  bcForLimits{sampleIndex} = bcValues(isfinite(bcValues));

  % Match the registered intragranular workflow: KAM is evaluated within
  % 2 degree reconstructed orientation domains, while values above 5
  % degrees are excluded from the neighbor average.
  [~, kamGrainId] = calcGrains(ebsdFull, "unitCell", ...
    "threshold", detectionDeg * degree);
  ebsdFull.grainId = kamGrainId;
  ebsdGrid = ebsdFull.gridify;
  kam = ebsdGrid.KAM("order", kamOrder, ...
    "threshold", kamCutoffDeg * degree);
  kamDeg = double(kam / degree);
  tiGrid = ebsdGrid("Ti-Hex");
  tiPhaseId = unique(double(tiGrid.phaseId));
  assert(isscalar(tiPhaseId));
  tiMask = reshape(double(ebsdGrid.phaseId), [], 1) == tiPhaseId;
  kamVector = reshape(kamDeg, [], 1);
  kamVector(~tiMask) = NaN;
  [gridX, gridY, gridIndices] = native_grid_indices(ebsdGrid);
  assert(isequal(gridX, xCoordinates{sampleIndex}) && ...
    isequal(gridY, yCoordinates{sampleIndex}), ...
    "Gridify changed registered map coordinates.");
  kamImages{sampleIndex} = scalar_image(kamVector, gridIndices, ...
    numel(gridY), numel(gridX));
  validKam = kamVector(isfinite(kamVector));
  assert(~isempty(validKam), "No finite Ti-Hex KAM values.");
  kamMean(sampleIndex) = mean(validKam);
  kamMedian(sampleIndex) = median(validKam);
  kamP90(sampleIndex) = prctile(validKam, 90);

  clear ebsdFull ebsdGrid tiEbsd tiGrid grains grainId kamGrainId kam
  clear innerBoundary outerBoundary boundaryA boundaryB boundaryC
end

allBc = vertcat(bcForLimits{:});
bcLimits = prctile(allBc, [1 99]);
if bcLimits(1) >= bcLimits(2)
  bcLimits = [min(allBc) max(allBc)];
end

pngPath = fullfile(outputDir, "vertical_ipf_gb_kam_matrix_raw.png");
tifPath = fullfile(outputDir, "vertical_ipf_gb_kam_matrix_raw.tif");
pdfPath = fullfile(outputDir, "vertical_ipf_gb_kam_matrix_raw.pdf");
render_matrix(ipfImages, bcImages, kamImages, xCoordinates, yCoordinates, ...
  boundary2To5, boundary5To15, boundary15Plus, catalog, ...
  tiCrystalSymmetry, bcLimits, kamDisplayDeg, kamMean, pngPath, ...
  tifPath, pdfPath);

totalBoundaryLength = boundaryLength2To5 + boundaryLength5To15 + ...
  boundaryLength15Plus;
sample = catalog.sample;
diameter_mm = catalog.diameter_mm;
cold_reduction_percent = catalog.cold_reduction_percent;
variant = catalog.variant;
input_path = catalog.input_path;
grain_detection_deg = repmat(detectionDeg, sampleCount, 1);
hagb_threshold_deg = repmat(hagbDeg, sampleCount, 1);
min_grain_pixels = repmat(minGrainPixels, sampleCount, 1);
kam_order = kamOrder * ones(sampleCount, 1);
kam_cutoff_deg = repmat(kamCutoffDeg, sampleCount, 1);
kam_mean_deg = kamMean;
kam_median_deg = kamMedian;
kam_p90_deg = kamP90;
lagb_2_5_length_fraction = boundaryLength2To5 ./ totalBoundaryLength;
lagb_5_15_length_fraction = boundaryLength5To15 ./ totalBoundaryLength;
hagb_ge15_length_fraction = boundaryLength15Plus ./ totalBoundaryLength;
summary = table(sample, diameter_mm, cold_reduction_percent, variant, ...
  input_path, grain_detection_deg, hagb_threshold_deg, min_grain_pixels, ...
  kam_order, kam_cutoff_deg, kam_mean_deg, kam_median_deg, kam_p90_deg, ...
  boundaryLength2To5, boundaryLength5To15, boundaryLength15Plus, ...
  lagb_2_5_length_fraction, lagb_5_15_length_fraction, ...
  hagb_ge15_length_fraction);
writetable(summary, fullfile(outputDir, ...
  "vertical_ipf_gb_kam_matrix_summary.csv"));
end

function render_matrix(ipfImages, bcImages, kamImages, xCoordinates, ...
    yCoordinates, boundary2To5, boundary5To15, boundary15Plus, ...
    catalog, tiCrystalSymmetry, bcLimits, kamLimits, kamMean, ...
    pngPath, tifPath, pdfPath)
figureHandle = figure("Visible", "off", "Color", "white", ...
  "Units", "inches", "Position", [0.25 0.25 8.1 16.2]);
cleanupFigure = onCleanup(@() close(figureHandle));

columnLeft = [0.12 0.405 0.69];
columnWidth = 0.27;
rowHeight = 0.135;
rowGap = 0.008;
topEdge = 0.955;
columnTitles = ["IPF-AD", "Grain-boundary classes", ...
  "KAM (order 1; 5 deg cutoff)"];
for columnIndex = 1:3
  annotation(figureHandle, "textbox", ...
    [columnLeft(columnIndex) 0.962 columnWidth 0.026], ...
    "String", columnTitles(columnIndex), "LineStyle", "none", ...
    "HorizontalAlignment", "center", "VerticalAlignment", "middle", ...
    "FontName", "Arial", "FontSize", 9, "FontWeight", "bold");
end

panelLetters = char('a' + (0:17));
for sampleIndex = 1:6
  panelBottom = topEdge - sampleIndex * rowHeight - ...
    (sampleIndex - 1) * rowGap;
  rowLabel = sprintf("%.2f%%\n%.2f mm", ...
    catalog.cold_reduction_percent(sampleIndex), ...
    catalog.diameter_mm(sampleIndex));
  annotation(figureHandle, "textbox", ...
    [0.006 panelBottom 0.105 rowHeight], "String", rowLabel, ...
    "LineStyle", "none", "HorizontalAlignment", "center", ...
    "VerticalAlignment", "middle", "FontName", "Arial", ...
    "FontSize", 8, "FontWeight", "bold");

  for columnIndex = 1:3
    axesHandle = axes(figureHandle, "Units", "normalized", ...
      "Position", [columnLeft(columnIndex), panelBottom, ...
      columnWidth, rowHeight]);
    switch columnIndex
      case 1
        imagesc(axesHandle, xCoordinates{sampleIndex}, ...
          yCoordinates{sampleIndex}, ipfImages{sampleIndex});
        format_map_axes(axesHandle);
        draw_boundary(axesHandle, boundary15Plus{sampleIndex, 1}, ...
          boundary15Plus{sampleIndex, 2}, [0.05 0.05 0.05], 0.35);
      case 2
        imagesc(axesHandle, xCoordinates{sampleIndex}, ...
          yCoordinates{sampleIndex}, bcImages{sampleIndex});
        format_map_axes(axesHandle);
        colormap(axesHandle, gray(256));
        clim(axesHandle, bcLimits);
        draw_boundary(axesHandle, boundary2To5{sampleIndex, 1}, ...
          boundary2To5{sampleIndex, 2}, [0.00 0.55 0.85], 0.45);
        draw_boundary(axesHandle, boundary5To15{sampleIndex, 1}, ...
          boundary5To15{sampleIndex, 2}, [0.95 0.55 0.00], 0.50);
        draw_boundary(axesHandle, boundary15Plus{sampleIndex, 1}, ...
          boundary15Plus{sampleIndex, 2}, [0.10 0.10 0.10], 0.55);
      case 3
        imageHandle = imagesc(axesHandle, xCoordinates{sampleIndex}, ...
          yCoordinates{sampleIndex}, kamImages{sampleIndex});
        set(imageHandle, "AlphaData", ...
          isfinite(kamImages{sampleIndex}));
        format_map_axes(axesHandle);
        set(axesHandle, "Color", [0.82 0.82 0.82]);
        colormap(axesHandle, turbo(256));
        clim(axesHandle, kamLimits);
        draw_boundary(axesHandle, boundary15Plus{sampleIndex, 1}, ...
          boundary15Plus{sampleIndex, 2}, [0.05 0.05 0.05], 0.30);
        text(axesHandle, 0.98, 0.04, ...
          sprintf("mean %.2f deg", kamMean(sampleIndex)), ...
          "Units", "normalized", "HorizontalAlignment", "right", ...
          "VerticalAlignment", "bottom", "FontName", "Arial", ...
          "FontSize", 6.5, "BackgroundColor", "white", ...
          "Margin", 1);
    end
    draw_scale_bar(axesHandle, xCoordinates{sampleIndex}, ...
      yCoordinates{sampleIndex}, 100);
    panelIndex = (sampleIndex - 1) * 3 + columnIndex;
    text(axesHandle, 0.015, 0.985, sprintf("(%c)", ...
      panelLetters(panelIndex)), "Units", "normalized", ...
      "HorizontalAlignment", "left", "VerticalAlignment", "top", ...
      "FontName", "Arial", "FontSize", 7, "FontWeight", "bold", ...
      "BackgroundColor", "white", "Margin", 1);
  end
end

draw_ipf_key(figureHandle, tiCrystalSymmetry, ...
  [columnLeft(1) 0.012 columnWidth 0.074]);
draw_boundary_legend(figureHandle, ...
  [columnLeft(2) 0.012 columnWidth 0.074]);
draw_kam_legend(figureHandle, ...
  [columnLeft(3) 0.020 columnWidth 0.050], kamLimits);
annotation(figureHandle, "textbox", [0.01 0.989 0.98 0.010], ...
  "String", "Raw EBSD | AD horizontal; TD/RD vertical", ...
  "LineStyle", "none", "HorizontalAlignment", "center", ...
  "VerticalAlignment", "top", "FontName", "Arial", ...
  "FontSize", 8, "FontWeight", "normal");

drawnow;
exportgraphics(figureHandle, pngPath, "Resolution", 600, ...
  "BackgroundColor", "white");
exportgraphics(figureHandle, tifPath, "Resolution", 600, ...
  "BackgroundColor", "white");
exportgraphics(figureHandle, pdfPath, "ContentType", "image", ...
  "Resolution", 600, "BackgroundColor", "white");
clear cleanupFigure
end

function imageData = build_ipf_image(ebsdFull, tiEbsd, nativeIndices, ...
    rowCount, columnCount)
imageData = repmat(reshape([0.82 0.82 0.82], 1, 1, 3), ...
  rowCount, columnCount, 1);
colorKey = ipfHSVKey(tiEbsd);
colorKey.inversePoleFigureDirection = xvector;
colors = colorKey.orientation2color(tiEbsd.orientations);
[found, fullRows] = ismember(double(tiEbsd.id(:)), double(ebsdFull.id(:)));
assert(all(found));
for channelIndex = 1:3
  channel = imageData(:, :, channelIndex);
  channel(nativeIndices(fullRows)) = colors(:, channelIndex);
  imageData(:, :, channelIndex) = channel;
end
end

function imageData = scalar_image(values, indices, rowCount, columnCount)
values = reshape(double(values), [], 1);
assert(numel(values) == numel(indices));
imageData = nan(rowCount, columnCount);
imageData(indices) = values;
end

function [xValues, yValues, linearIndices] = native_grid_indices(ebsd)
xCoordinates = double(ebsd.x(:));
yCoordinates = double(ebsd.y(:));
xValues = unique(xCoordinates);
yValues = unique(yCoordinates);
[xFound, xIndex] = ismember(xCoordinates, xValues);
[yFound, yIndex] = ismember(yCoordinates, yValues);
assert(all(xFound & yFound));
linearIndices = sub2ind([numel(yValues), numel(xValues)], yIndex, xIndex);
assert(numel(unique(linearIndices)) == length(ebsd));
end

function [xCoordinates, yCoordinates] = boundary_coordinates(boundary)
if isempty(boundary)
  xCoordinates = nan;
  yCoordinates = nan;
  return
end
vertices = boundary.allV.xyz;
faces = boundary.F;
xCoordinates = reshape(vertices(faces.', 1), 2, []);
yCoordinates = reshape(vertices(faces.', 2), 2, []);
separator = nan(1, size(xCoordinates, 2));
xCoordinates = reshape([xCoordinates; separator], [], 1);
yCoordinates = reshape([yCoordinates; separator], [], 1);
end

function lengthUm = boundary_length(boundary)
if isempty(boundary)
  lengthUm = 0;
else
  lengthUm = sum(double(boundary.segLength));
end
end

function format_map_axes(axesHandle)
axis(axesHandle, "image");
set(axesHandle, "YDir", "normal", "XTick", [], "YTick", [], ...
  "Box", "on", "LineWidth", 0.45, "FontName", "Arial");
hold(axesHandle, "on");
end

function draw_boundary(axesHandle, xCoordinates, yCoordinates, color, width)
line(axesHandle, xCoordinates, yCoordinates, "Color", color, ...
  "LineWidth", width);
end

function draw_scale_bar(axesHandle, xValues, yValues, lengthUm)
xRange = max(xValues) - min(xValues);
yRange = max(yValues) - min(yValues);
xStart = min(xValues) + 0.055 * xRange;
yPosition = min(yValues) + 0.065 * yRange;
line(axesHandle, [xStart xStart + lengthUm], [yPosition yPosition], ...
  "Color", "black", "LineWidth", 1.6);
text(axesHandle, xStart + lengthUm / 2, yPosition + 0.018 * yRange, ...
  sprintf("%d um", lengthUm), "HorizontalAlignment", "center", ...
  "VerticalAlignment", "bottom", "FontName", "Arial", ...
  "FontSize", 5.8, "BackgroundColor", "white", "Margin", 0.5);
end

function draw_ipf_key(figureHandle, crystalSymmetry, position)
axesHandle = axes(figureHandle, "Units", "normalized", ...
  "Position", position);
colorKey = ipfHSVKey(crystalSymmetry);
colorKey.inversePoleFigureDirection = xvector;
plot(colorKey, "parent", axesHandle, "noTitle");
title(axesHandle, "Ti-Hex IPF key || AD", "FontName", "Arial", ...
  "FontSize", 7, "FontWeight", "normal");
end

function draw_boundary_legend(figureHandle, position)
axesHandle = axes(figureHandle, "Units", "normalized", ...
  "Position", position);
axis(axesHandle, [0 1 0 1]);
axis(axesHandle, "off");
colors = [0.00 0.55 0.85; 0.95 0.55 0.00; 0.10 0.10 0.10];
labels = ["2-<5 deg", "5-<15 deg", ">=15 deg"];
for index = 1:3
  y = 0.82 - (index - 1) * 0.31;
  line(axesHandle, [0.06 0.30], [y y], "Color", colors(index, :), ...
    "LineWidth", 1.5);
  text(axesHandle, 0.36, y, labels(index), "FontName", "Arial", ...
    "FontSize", 7, "VerticalAlignment", "middle");
end
title(axesHandle, "Ti-Hex boundary classes", "FontName", "Arial", ...
  "FontSize", 7, "FontWeight", "normal");
end

function draw_kam_legend(figureHandle, position, limits)
axesHandle = axes(figureHandle, "Units", "normalized", ...
  "Position", position);
imagesc(axesHandle, linspace(limits(1), limits(2), 256), 1, 1:256);
set(axesHandle, "YTick", [], "XTick", limits, "XLim", limits, ...
  "FontName", "Arial", "FontSize", 7, "Box", "on");
colormap(axesHandle, turbo(256));
xlabel(axesHandle, "KAM (deg)", "FontName", "Arial", "FontSize", 7);
end
