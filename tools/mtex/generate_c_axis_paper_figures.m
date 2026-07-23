function metadata = generate_c_axis_paper_figures(scanRoot,outputRoot)
arguments
  scanRoot (1,1) string
  outputRoot (1,1) string
end
assert(isfolder(scanRoot));
assert(~isempty(which("EBSD")));
textureDir = fullfile(outputRoot,"05_texture");
summaryPath = fullfile(textureDir,"texture_summary.csv");
assert(isfile(summaryPath));
summary = readtable(summaryPath,"TextType","string", ...
  "VariableNamingRule","preserve");
prepared = prepare_c_axis_figure_rows(summary);
catalog = comprehensive_ebsd_catalog(scanRoot);
assert(height(catalog) == 12);
parameters = comprehensive_ebsd_output_contract().parameters;
kernel = SO3DeLaValleePoussinKernel("halfwidth", ...
  parameters.texture_kernel_halfwidth_deg * degree);
odfs = cell(12,1);
sharedMaximum = 0;
for scanIndex = 1:12
  [ebsdFull,~] = load_comprehensive_ebsd_scan(catalog(scanIndex,:));
  tiEbsd = ebsdFull("Ti-Hex");
  rbfOdf = calcDensity(tiEbsd.orientations,"kernel",kernel, ...
    "weights",ones(length(tiEbsd),1),"silent");
  rbfOdf = rbfOdf / double(mean(rbfOdf));
  odfs{scanIndex} = SO3FunHarmonic(rbfOdf, ...
    "bandwidth",kernel.bandwidth);
  cAxis = Miller(0,0,0,1,tiEbsd.CS);
  density = calcPDF(odfs{scanIndex},cAxis,[],"antipodal");
  sharedMaximum = max(sharedMaximum,double(max(density, ...
    "resolution",parameters.texture_grid_resolution_deg*degree)));
  clear ebsdFull tiEbsd rbfOdf density
end
assert(isfinite(sharedMaximum) && sharedMaximum > 0);

previousAnnotations = getMTEXpref("pfAnnotations");
restoreAnnotations = onCleanup(@() setMTEXpref( ...
  "pfAnnotations",previousAnnotations));
pfAnnotations = @(varargin) text( ...
  [vector3d.X,vector3d.Y,vector3d.Z], ...
  {"AD","TD-RD","ND"}, ...
  "BackgroundColor","w","tag","axesLabels",varargin{:});
setMTEXpref("pfAnnotations",pfAnnotations);

render_pole_montage(odfs,catalog,"raw",sharedMaximum, ...
  fullfile(textureDir,"c_axis_pole_figures_raw.png"),parameters);
render_pole_montage(odfs,catalog,"denoised",sharedMaximum, ...
  fullfile(textureDir,"c_axis_pole_figures_denoised.png"),parameters);
rawPlot = render_mean_trend( ...
  prepared.raw,"raw",prepared.common_y_limits_deg, ...
  fullfile(textureDir,"c_axis_mean_orientation_raw.png"));
denoisedPlot = render_mean_trend(prepared.denoised,"denoised", ...
  prepared.common_y_limits_deg, ...
  fullfile(textureDir,"c_axis_mean_orientation_denoised.png"));

metadata = build_c_axis_figure_metadata(sharedMaximum, ...
  prepared.common_y_limits_deg,rawPlot,denoisedPlot);
clear restoreAnnotations
end

function render_pole_montage(odfs,catalog,variantName, ...
  colorMaximum,outputPath,parameters)
rows = find(catalog.variant == variantName);
assert(numel(rows) == 6);
[~,order] = sort(catalog.cold_reduction_percent(rows));
rows = rows(order);
figureHandle = figure("Visible","off","Color","white", ...
  "Position",[50 50 1800 1180]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle,2,3,"Padding","compact", ...
  "TileSpacing","compact");
for panelIndex = 1:6
  axesHandle = nexttile(layout);
  odf = odfs{rows(panelIndex)};
  cAxis = Miller(0,0,0,1,odf.CS);
  plotPDF(odf,cAxis,"antipodal","earea","contourf","silent", ...
    "parent",axesHandle, ...
    "resolution",parameters.texture_grid_resolution_deg*degree, ...
    "colorRange",[0 colorMaximum]);
  title(axesHandle,sprintf("%.2f%%", ...
    catalog.cold_reduction_percent(rows(panelIndex))), ...
    "Interpreter","none","FontWeight","normal");
end
colorbarHandle = colorbar(axesHandle,"eastoutside");
colorbarHandle.Label.String = "MRD";
title(layout,sprintf( ...
  "Alpha-Ti {0001} c-axis pole figures | %s | AD horizontal", ...
  variantName),"Interpreter","none");
exportgraphics(figureHandle,outputPath,"Resolution",300, ...
  "BackgroundColor","white");
clear cleanupFigure
end

function plotData = render_mean_trend( ...
  rows,variantName,yLimits,outputPath)
figureHandle = figure("Visible","off","Color","white", ...
  "Position",[100 100 1200 760]);
cleanupFigure = onCleanup(@() close(figureHandle));
axesHandle = axes(figureHandle);
hold(axesHandle,"on");
x = rows.cold_reduction_percent;
meanDeg = rows.c_axis_ad_mean_deg;
p10Deg = rows.c_axis_ad_p10_deg;
p90Deg = rows.c_axis_ad_p90_deg;
fill(axesHandle,[x;flipud(x)], ...
  [p10Deg;flipud(p90Deg)], ...
  [0.75 0.84 0.94],"EdgeColor","none","FaceAlpha",0.45, ...
  "DisplayName","P10-P90");
plot(axesHandle,x,meanDeg,"-o", ...
  "LineWidth",1.8,"MarkerSize",7, ...
  "DisplayName","Mean c-axis-AD angle");
for pointIndex = 1:height(rows)
  text(axesHandle,x(pointIndex),meanDeg(pointIndex), ...
    sprintf("  %.2f",meanDeg(pointIndex)), ...
    "VerticalAlignment","bottom");
end
xlim(axesHandle,[0 50]);
ylim(axesHandle,yLimits);
xlabel(axesHandle,"Cold reduction (%)");
ylabel(axesHandle,"c-axis-AD acute angle (deg)");
title(axesHandle,sprintf( ...
  "Mean alpha-Ti c-axis orientation | %s | AD horizontal", ...
  variantName),"Interpreter","none","FontWeight","normal");
grid(axesHandle,"on");
legend(axesHandle,"Location","best");
exportgraphics(figureHandle,outputPath,"Resolution",300, ...
  "BackgroundColor","white");
plotData = struct();
plotData.x = x;
plotData.mean_deg = meanDeg;
plotData.p10_deg = p10Deg;
plotData.p90_deg = p90Deg;
clear cleanupFigure
end
