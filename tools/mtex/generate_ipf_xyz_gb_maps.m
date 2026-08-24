function summary = generate_ipf_xyz_gb_maps(scanRoot,outputDir)
%GENERATE_IPF_XYZ_GB_MAPS Export raw IPF-X/Y/Z maps with Ti-Hex boundaries.

arguments
  scanRoot (1,1) string
  outputDir (1,1) string
end

assert(isfolder(scanRoot),"EBSD scan folder not found: %s",scanRoot);
assert(~isempty(which("EBSD")),"MTEX must be loaded.");
if ~isfolder(outputDir)
  mkdir(outputDir);
end

catalog = comprehensive_ebsd_catalog(scanRoot);
catalog = catalog(catalog.variant == "raw",:);
[~,order] = sort(catalog.cold_reduction_percent);
catalog = catalog(order,:);
assert(height(catalog) == 6);

detectionDeg = 2;
classificationDeg = 10;
directions = {xvector,yvector,zvector};
directionTitles = ["IPF-X / AD","IPF-Y / TD-RD","IPF-Z / ND"];
directionSuffixes = ["x_ad","y_td_rd","z_nd"];
sampleCount = height(catalog);
directionCount = numel(directions);
mapData = cell(sampleCount,directionCount);
xCoordinates = cell(sampleCount,1);
yCoordinates = cell(sampleCount,1);
lowBoundaryX = cell(sampleCount,1);
lowBoundaryY = cell(sampleCount,1);
highBoundaryX = cell(sampleCount,1);
highBoundaryY = cell(sampleCount,1);
pointCount = zeros(sampleCount,1);
lowLength = zeros(sampleCount,1);
highLength = zeros(sampleCount,1);
tiCrystalSymmetry = [];

for sampleIndex = 1:sampleCount
  fprintf("IPF_XYZ_GB sample=%s\n",catalog.sample(sampleIndex));
  [ebsdFull,~] = load_comprehensive_ebsd_scan(catalog(sampleIndex,:));
  tiEbsd = ebsdFull("Ti-Hex");
  assert(~isempty(tiEbsd));
  if isempty(tiCrystalSymmetry)
    tiCrystalSymmetry = tiEbsd.CS;
  end
  pointCount(sampleIndex) = length(tiEbsd);

  [grains,grainId] = calcGrains(ebsdFull,"unitCell", ...
    "threshold",[classificationDeg detectionDeg] * degree);
  ebsdFull.grainId = grainId;
  lowSource = grains.innerBoundary("Ti-Hex","Ti-Hex");
  highSource = grains.boundary("Ti-Hex","Ti-Hex");
  lowAngles = double(angle(lowSource.misorientation) / degree);
  highAngles = double(angle(highSource.misorientation) / degree);
  lowBoundary = lowSource(lowAngles >= detectionDeg & ...
    lowAngles < classificationDeg);
  highBoundary = highSource(highAngles >= classificationDeg);
  assert(~isempty(lowBoundary) && ~isempty(highBoundary));
  lowLength(sampleIndex) = sum(double(lowBoundary.segLength));
  highLength(sampleIndex) = sum(double(highBoundary.segLength));
  [lowBoundaryX{sampleIndex},lowBoundaryY{sampleIndex}] = ...
    boundary_coordinates(lowBoundary);
  [highBoundaryX{sampleIndex},highBoundaryY{sampleIndex}] = ...
    boundary_coordinates(highBoundary);

  [xCoordinates{sampleIndex},yCoordinates{sampleIndex},gridIndices] = ...
    native_grid_indices(ebsdFull);
  for directionIndex = 1:directionCount
    mapData{sampleIndex,directionIndex} = build_ipf_image( ...
      ebsdFull,directions{directionIndex},gridIndices, ...
      numel(yCoordinates{sampleIndex}),numel(xCoordinates{sampleIndex}));
    directionPath = fullfile(outputDir,catalog.sample(sampleIndex) + ...
      "_ipf_" + directionSuffixes(directionIndex) + "_gb.png");
    render_standalone_direction( ...
      mapData{sampleIndex,directionIndex}, ...
      xCoordinates{sampleIndex},yCoordinates{sampleIndex}, ...
      lowBoundaryX{sampleIndex},lowBoundaryY{sampleIndex}, ...
      highBoundaryX{sampleIndex},highBoundaryY{sampleIndex}, ...
      tiEbsd,directions{directionIndex}, ...
      directionTitles(directionIndex),catalog(sampleIndex,:), ...
      lowLength(sampleIndex),highLength(sampleIndex),directionPath);
  end

  samplePath = fullfile(outputDir,catalog.sample(sampleIndex) + ...
    "_ipf_xyz_gb.png");
  render_sample_figure(mapData(sampleIndex,:), ...
    xCoordinates{sampleIndex},yCoordinates{sampleIndex}, ...
    lowBoundaryX{sampleIndex},lowBoundaryY{sampleIndex}, ...
    highBoundaryX{sampleIndex},highBoundaryY{sampleIndex}, ...
    catalog(sampleIndex,:),lowLength(sampleIndex), ...
    highLength(sampleIndex),directionTitles,samplePath);
  clear ebsdFull tiEbsd grains grainId lowSource highSource
  clear lowBoundary highBoundary
end

function render_standalone_direction(imageData,xValues,yValues, ...
    lowX,lowY,highX,highY,tiEbsd,direction,directionTitle, ...
    catalogRow,lowLength,highLength,outputPath)
figureHandle = figure("Visible","off","Color","white", ...
  "Units","inches","Position",[0.25 0.25 8.5 5]);
cleanupFigure = onCleanup(@() close(figureHandle));

mapAxes = axes(figureHandle,"Units","normalized", ...
  "Position",[0.04 0.08 0.52 0.82]);
draw_map_panel(mapAxes,imageData,xValues,yValues, ...
  lowX,lowY,highX,highY);
title(mapAxes,directionTitle,"Interpreter","none", ...
  "FontWeight","bold","FontSize",11);

keyAxes = axes(figureHandle,"Units","normalized", ...
  "Position",[0.64 0.56 0.27 0.30]);
colorKey = ipfHSVKey(tiEbsd);
colorKey.inversePoleFigureDirection = direction;
plot(colorKey,"parent",keyAxes,"noTitle");
title(keyAxes,"Ti-Hex IPF key","Interpreter","none", ...
  "FontWeight","normal","FontSize",9);

boundaryAxes = axes(figureHandle,"Units","normalized", ...
  "Position",[0.63 0.31 0.31 0.17]);
draw_boundary_legend(boundaryAxes,lowLength,highLength);

coordinateAxes = axes(figureHandle,"Units","normalized", ...
  "Position",[0.63 0.06 0.31 0.18]);
draw_coordinate_legend(coordinateAxes);

sgtitle(figureHandle,sprintf("%.2f mm | %.2f%% cold reduction | raw", ...
  catalogRow.diameter_mm,catalogRow.cold_reduction_percent), ...
  "Interpreter","none","FontWeight","bold","FontSize",12);
drawnow;
exportgraphics(figureHandle,outputPath,"Resolution",600, ...
  "BackgroundColor","white");
clear cleanupFigure
end

function draw_boundary_legend(axesHandle,lowLength,highLength)
totalLength = lowLength + highLength;
axis(axesHandle,[0 1 0 1]);
axis(axesHandle,"off");
text(axesHandle,0,0.95,"Ti-Hex grain boundaries", ...
  "FontWeight","bold","FontSize",8,"VerticalAlignment","top");
line(axesHandle,[0.04 0.28],[0.60 0.60], ...
  "Color",[0.55 0.55 0.55],"LineWidth",1.0);
text(axesHandle,0.34,0.60,sprintf("2-<10 deg  %.1f%%", ...
  100 * lowLength / totalLength),"FontSize",8, ...
  "VerticalAlignment","middle");
line(axesHandle,[0.04 0.28],[0.26 0.26], ...
  "Color","black","LineWidth",1.8);
text(axesHandle,0.34,0.26,sprintf(">=10 deg  %.1f%%", ...
  100 * highLength / totalLength),"FontSize",8, ...
  "VerticalAlignment","middle");
end

function draw_coordinate_legend(axesHandle)
axis(axesHandle,[0 1 0 1]);
axis(axesHandle,"off");
hold(axesHandle,"on");
origin = [0.18 0.22];
quiver(axesHandle,origin(1),origin(2),0.48,0,0, ...
  "Color","black","LineWidth",1.2,"MaxHeadSize",0.25);
quiver(axesHandle,origin(1),origin(2),0,0.50,0, ...
  "Color","black","LineWidth",1.2,"MaxHeadSize",0.25);
text(axesHandle,0.71,0.22,"X = AD (axial)", ...
  "FontSize",8,"VerticalAlignment","middle");
text(axesHandle,0.18,0.78,"Y = TD/RD", ...
  "FontSize",8,"HorizontalAlignment","center");
text(axesHandle,0.18,0.04,"Z = ND (out of plane)", ...
  "FontSize",8,"HorizontalAlignment","center");
end

render_matrix_figure(mapData,xCoordinates,yCoordinates, ...
  lowBoundaryX,lowBoundaryY,highBoundaryX,highBoundaryY,catalog, ...
  directionTitles,fullfile(outputDir,"ipf_xyz_gb_raw_matrix.png"), ...
  fullfile(outputDir,"ipf_xyz_gb_raw_matrix.pdf"));
render_axial_montage(mapData(:,1),xCoordinates,yCoordinates, ...
  lowBoundaryX,lowBoundaryY,highBoundaryX,highBoundaryY,catalog, ...
  lowLength,highLength,tiCrystalSymmetry, ...
  fullfile(outputDir,"ipf_x_ad_six_state_montage.png"), ...
  fullfile(outputDir,"ipf_x_ad_six_state_montage.pdf"));

totalLength = lowLength + highLength;
sample = catalog.sample;
diameter_mm = catalog.diameter_mm;
cold_reduction_percent = catalog.cold_reduction_percent;
variant = catalog.variant;
input_path = catalog.input_path;
valid_ti_hex_point_count = pointCount;
grain_detection_deg = repmat(detectionDeg,sampleCount,1);
boundary_classification_deg = repmat(classificationDeg,sampleCount,1);
lagb_2_10_length_um = lowLength;
hagb_ge10_length_um = highLength;
lagb_2_10_length_fraction = lowLength ./ totalLength;
hagb_ge10_length_fraction = highLength ./ totalLength;
coordinate_x = repmat("AD",sampleCount,1);
coordinate_y = repmat("TD/RD",sampleCount,1);
coordinate_z = repmat("ND",sampleCount,1);
summary = table(sample,diameter_mm,cold_reduction_percent,variant, ...
  input_path,valid_ti_hex_point_count,grain_detection_deg, ...
  boundary_classification_deg,lagb_2_10_length_um, ...
  hagb_ge10_length_um,lagb_2_10_length_fraction, ...
  hagb_ge10_length_fraction,coordinate_x,coordinate_y,coordinate_z);
writetable(summary,fullfile(outputDir,"ipf_xyz_gb_summary.csv"));
end

function render_axial_montage(images,xValues,yValues,lowX,lowY, ...
    highX,highY,catalog,lowLength,highLength,tiCrystalSymmetry, ...
    pngPath,pdfPath)
figureHandle = figure("Visible","off","Color","white", ...
  "Units","inches","Position",[0.25 0.25 13 8]);
cleanupFigure = onCleanup(@() close(figureHandle));
panelLabels = ["(a)","(b)","(c)","(d)","(e)","(f)"];
panelLeft = [0.02 0.28 0.54];
panelBottom = [0.53 0.07];
panelWidth = 0.235;
panelHeight = 0.37;

for sampleIndex = 1:6
  rowIndex = 1 + (sampleIndex > 3);
  columnIndex = mod(sampleIndex - 1,3) + 1;
  axesHandle = axes(figureHandle,"Units","normalized", ...
    "Position",[panelLeft(columnIndex),panelBottom(rowIndex), ...
    panelWidth,panelHeight]);
  draw_map_panel(axesHandle,images{sampleIndex}, ...
    xValues{sampleIndex},yValues{sampleIndex}, ...
    lowX{sampleIndex},lowY{sampleIndex}, ...
    highX{sampleIndex},highY{sampleIndex});
  totalLength = lowLength(sampleIndex) + highLength(sampleIndex);
  title(axesHandle,sprintf( ...
    "%s %.2f mm | %.2f%%\n2-<10 deg %.1f%% | >=10 deg %.1f%%", ...
    panelLabels(sampleIndex),catalog.diameter_mm(sampleIndex), ...
    catalog.cold_reduction_percent(sampleIndex), ...
    100 * lowLength(sampleIndex) / totalLength, ...
    100 * highLength(sampleIndex) / totalLength), ...
    "Interpreter","none","FontWeight","normal","FontSize",8);
end

keyAxes = axes(figureHandle,"Units","normalized", ...
  "Position",[0.81 0.60 0.16 0.24]);
colorKey = ipfHSVKey(tiCrystalSymmetry);
colorKey.inversePoleFigureDirection = xvector;
plot(colorKey,"parent",keyAxes,"noTitle");
title(keyAxes,"Ti-Hex IPF key","Interpreter","none", ...
  "FontWeight","normal","FontSize",9);

boundaryAxes = axes(figureHandle,"Units","normalized", ...
  "Position",[0.81 0.36 0.17 0.16]);
draw_shared_boundary_legend(boundaryAxes);

coordinateAxes = axes(figureHandle,"Units","normalized", ...
  "Position",[0.80 0.08 0.18 0.20]);
draw_axial_coordinate_legend(coordinateAxes);

sgtitle(figureHandle, ...
  "Alpha-Ti IPF-X maps (coloring || AD, axial direction)", ...
  "Interpreter","none","FontWeight","bold","FontSize",12);
drawnow;
exportgraphics(figureHandle,pngPath,"Resolution",600, ...
  "BackgroundColor","white");
exportgraphics(figureHandle,pdfPath,"ContentType","image", ...
  "Resolution",600,"BackgroundColor","white");
clear cleanupFigure
end

function draw_shared_boundary_legend(axesHandle)
axis(axesHandle,[0 1 0 1]);
axis(axesHandle,"off");
text(axesHandle,0,0.95,"Ti-Hex grain boundaries", ...
  "FontWeight","bold","FontSize",8,"VerticalAlignment","top");
line(axesHandle,[0.04 0.32],[0.57 0.57], ...
  "Color",[0.55 0.55 0.55],"LineWidth",1.0);
text(axesHandle,0.38,0.57,"2-<10 deg", ...
  "FontSize",8,"VerticalAlignment","middle");
line(axesHandle,[0.04 0.32],[0.24 0.24], ...
  "Color","black","LineWidth",1.8);
text(axesHandle,0.38,0.24,">=10 deg", ...
  "FontSize",8,"VerticalAlignment","middle");
end

function draw_axial_coordinate_legend(axesHandle)
axis(axesHandle,[0 1 0 1]);
axis(axesHandle,"off");
hold(axesHandle,"on");
origin = [0.18 0.22];
quiver(axesHandle,origin(1),origin(2),0.48,0,0, ...
  "Color","black","LineWidth",1.2,"MaxHeadSize",0.25);
quiver(axesHandle,origin(1),origin(2),0,0.50,0, ...
  "Color","black","LineWidth",1.2,"MaxHeadSize",0.25);
text(axesHandle,0.71,0.22,"X = AD (axial)", ...
  "FontSize",8,"VerticalAlignment","middle");
text(axesHandle,0.18,0.78,"Y = TD/RD", ...
  "FontSize",8,"HorizontalAlignment","center");
text(axesHandle,0.18,0.04,"Z = ND (out of plane)", ...
  "FontSize",8,"HorizontalAlignment","center");
end

function imageData = build_ipf_image(ebsdFull,direction,gridIndices, ...
    rowCount,columnCount)
imageData = uint8(255 * ones(rowCount,columnCount,3));
phaseNames = ["Ti-Hex","Titanium cubic"];
for phaseIndex = 1:numel(phaseNames)
  phaseEbsd = ebsdFull(phaseNames(phaseIndex));
  if isempty(phaseEbsd)
    continue
  end
  colorKey = ipfHSVKey(phaseEbsd);
  colorKey.inversePoleFigureDirection = direction;
  phaseColors = colorKey.orientation2color(phaseEbsd.orientations);
  [found,fullRows] = ismember(double(phaseEbsd.id(:)), ...
    double(ebsdFull.id(:)));
  assert(all(found));
  for colorIndex = 1:3
    channel = imageData(:,:,colorIndex);
    channel(gridIndices(fullRows)) = uint8(round( ...
      255 * phaseColors(:,colorIndex)));
    imageData(:,:,colorIndex) = channel;
  end
end
end

function render_sample_figure(images,xValues,yValues,lowX,lowY, ...
    highX,highY,catalogRow,lowLength,highLength,directionTitles, ...
    outputPath)
figureHandle = figure("Visible","off","Color","white", ...
  "Units","inches","Position",[0.25 0.25 12 4.25]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle,1,3,"Padding","compact", ...
  "TileSpacing","compact");
for directionIndex = 1:3
  axesHandle = nexttile(layout);
  draw_map_panel(axesHandle,images{directionIndex},xValues,yValues, ...
    lowX,lowY,highX,highY);
  title(axesHandle,directionTitles(directionIndex), ...
    "Interpreter","none","FontWeight","normal");
end
totalLength = lowLength + highLength;
title(layout,sprintf("%.2f mm | %.2f%% cold reduction | raw", ...
  catalogRow.diameter_mm,catalogRow.cold_reduction_percent), ...
  "Interpreter","none","FontWeight","bold");
xlabel(layout,sprintf( ...
  "Ti-Hex boundaries by length: 2-<10 deg %.1f%% | >=10 deg %.1f%%", ...
  100 * lowLength / totalLength,100 * highLength / totalLength));
exportgraphics(figureHandle,outputPath,"Resolution",600, ...
  "BackgroundColor","white");
clear cleanupFigure
end

function render_matrix_figure(images,xValues,yValues,lowX,lowY, ...
    highX,highY,catalog,directionTitles,pngPath,pdfPath)
figureHandle = figure("Visible","off","Color","white", ...
  "Units","inches","Position",[0.25 0.25 10 16]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle,6,3,"Padding","compact", ...
  "TileSpacing","compact");
for sampleIndex = 1:6
  for directionIndex = 1:3
    axesHandle = nexttile(layout);
    draw_map_panel(axesHandle,images{sampleIndex,directionIndex}, ...
      xValues{sampleIndex},yValues{sampleIndex}, ...
      lowX{sampleIndex},lowY{sampleIndex}, ...
      highX{sampleIndex},highY{sampleIndex});
    title(axesHandle,sprintf("%.2f mm | %.2f%% | %s", ...
      catalog.diameter_mm(sampleIndex), ...
      catalog.cold_reduction_percent(sampleIndex), ...
      directionTitles(directionIndex)), ...
      "Interpreter","none","FontWeight","normal","FontSize",8);
  end
end
title(layout,"Raw EBSD IPF maps with Ti-Hex grain boundaries", ...
  "Interpreter","none","FontWeight","bold");
xlabel(layout, ...
  "Gray: 2-<10 deg | Black: >=10 deg | X=AD, Y=TD/RD, Z=ND");
exportgraphics(figureHandle,pngPath,"Resolution",600, ...
  "BackgroundColor","white");
exportgraphics(figureHandle,pdfPath,"ContentType","image", ...
  "Resolution",600,"BackgroundColor","white");
clear cleanupFigure
end

function draw_map_panel(axesHandle,imageData,xValues,yValues, ...
    lowX,lowY,highX,highY)
imagesc(axesHandle,xValues,yValues,imageData);
axis(axesHandle,"image");
set(axesHandle,"YDir","normal","XTick",[],"YTick",[], ...
  "Box","on");
hold(axesHandle,"on");
line(axesHandle,lowX,lowY,"Color",[0.55 0.55 0.55], ...
  "LineWidth",0.25);
line(axesHandle,highX,highY,"Color","black","LineWidth",0.55);
draw_scale_bar(axesHandle,xValues,yValues,100);
end

function draw_scale_bar(axesHandle,xValues,yValues,lengthUm)
xRange = max(xValues) - min(xValues);
yRange = max(yValues) - min(yValues);
xStart = min(xValues) + 0.06 * xRange;
yPosition = min(yValues) + 0.07 * yRange;
line(axesHandle,[xStart xStart + lengthUm],[yPosition yPosition], ...
  "Color","black","LineWidth",2.0);
text(axesHandle,xStart + lengthUm / 2,yPosition + 0.025 * yRange, ...
  sprintf("%d um",lengthUm),"HorizontalAlignment","center", ...
  "VerticalAlignment","bottom","FontSize",7, ...
  "BackgroundColor","white");
end

function [xValues,yValues,linearIndices] = native_grid_indices(ebsdFull)
xCoordinates = double(ebsdFull.x(:));
yCoordinates = double(ebsdFull.y(:));
xValues = unique(xCoordinates);
yValues = unique(yCoordinates);
[xFound,xIndex] = ismember(xCoordinates,xValues);
[yFound,yIndex] = ismember(yCoordinates,yValues);
assert(all(xFound & yFound));
linearIndices = sub2ind([numel(yValues),numel(xValues)], ...
  yIndex,xIndex);
assert(numel(unique(linearIndices)) == length(ebsdFull));
end

function [xCoordinates,yCoordinates] = boundary_coordinates(boundary)
vertices = boundary.allV.xyz;
faces = boundary.F;
xCoordinates = reshape(vertices(faces.',1),2,[]);
yCoordinates = reshape(vertices(faces.',2),2,[]);
separator = nan(1,size(xCoordinates,2));
xCoordinates = reshape([xCoordinates;separator],[],1);
yCoordinates = reshape([yCoordinates;separator],[],1);
end
