function summary = generate_c_axis_pole_figures(scanRoot, outputDir)
%GENERATE_C_AXIS_POLE_FIGURES Create area-weighted alpha-Ti {0001} pole figures.

arguments
  scanRoot (1,1) string
  outputDir (1,1) string
end

assert(isfolder(scanRoot), "EBSD scan folder not found: %s", scanRoot);
assert(~isempty(which("EBSD")), "MTEX is not loaded in the current MATLAB session.");

if ~isfolder(outputDir)
  mkdir(outputDir);
end

folderNames = ["d7"; "d6_48"; "d6_02"; "d5_6"; "d5_25"; "d5"];
inputNames = [
  "ebsd_sample_7_map_15.ctf";
  "ebsd_sample_648_map_13.ctf";
  "ebsd_sample_602_map_11.ctf";
  "ebsd_sample_56_map_9.ctf";
  "ebsd_sample_525_map_7.ctf";
  "ebsd_sample_5_map_3.ctf"
];
sampleNames = ["7d"; "6.48d"; "6.02d"; "5.6d"; "5.25d"; "5d"];
coldReductionPercent = [0; 14.31; 26.04; 36.00; 43.75; 48.98];

n = numel(sampleNames);
odfs = cell(n,1);
cAxes = cell(n,1);
totalPixelCount = zeros(n,1);
alphaPixelCount = zeros(n,1);
maxMRD = zeros(n,1);
peakDirection = cell(n,1);

kernel = SO3DeLaValleePoussinKernel("halfwidth", 5 * degree);

for i = 1:n
  inputFile = fullfile(scanRoot, folderNames(i), inputNames(i));
  assert(isfile(inputFile), "CTF file not found: %s", inputFile);

  fprintf("IMPORT=%s\n", inputFile);
  ebsd = EBSD.load(inputFile, "convertEuler2SpatialReferenceFrame");
  alpha = ebsd("Ti-Hex");

  assert(~isempty(alpha), "Ti-Hex phase missing in %s", inputFile);
  totalPixelCount(i) = length(ebsd);
  alphaPixelCount(i) = length(alpha);

  fprintf("ODF=%s ALPHA_PIXELS=%d\n", sampleNames(i), alphaPixelCount(i));
  odfs{i} = calcDensity(alpha.orientations, "kernel", kernel);
  cAxes{i} = Miller(0, 0, 0, 1, alpha.CS);
  pdf = odfs{i}.calcPDF(cAxes{i}, [], "antipodal");
  [maxMRD(i), peakDirection{i}] = max(pdf);

  clear ebsd alpha pdf
end

colorMaximum = ceil(max(maxMRD));
if colorMaximum <= 0
  colorMaximum = 1;
end

previousVisibility = get(groot, "defaultFigureVisible");
restoreVisibility = onCleanup(@() set(groot, ...
  "defaultFigureVisible", previousVisibility));
set(groot, "defaultFigureVisible", "off");

pfAnnotations = @(varargin) text( ...
  [vector3d.X, vector3d.Y, vector3d.Z], {"AD", "RD", "ND"}, ...
  "BackgroundColor", "w", "tag", "axesLabels", varargin{:});
setMTEXpref("pfAnnotations", pfAnnotations);

outputFiles = strings(n,1);
peakX = zeros(n,1);
peakY = zeros(n,1);
peakZ = zeros(n,1);
cAxisToADAngleDeg = zeros(n,1);

for i = 1:n
  close all
  plotPDF(odfs{i}, cAxes{i}, "antipodal", "contourf");
  setColorRange([0, colorMaximum], "current");
  mtexColorbar("title", "m.r.d.");
  mtexTitle(sampleNames(i));

  outputFiles(i) = fullfile(outputDir, sampleNames(i) + ".png");
  figureHandle = gcf;
  set(figureHandle, "Position", [100, 100, 900, 900]);
  axesHandles = findall(figureHandle, "Type", "axes");
  for ax = reshape(axesHandles, 1, [])
    if isprop(ax, "Toolbar")
      ax.Toolbar = [];
    end
  end
  drawnow;
  exportgraphics(figureHandle, char(outputFiles(i)), ...
    "Resolution", 300, "BackgroundColor", "white");
  renderedImage = imread(outputFiles(i));
  imwrite(renderedImage, outputFiles(i), "png");
  assert(isfile(outputFiles(i)), "Output was not created: %s", outputFiles(i));

  peakX(i) = peakDirection{i}.x;
  peakY(i) = peakDirection{i}.y;
  peakZ(i) = peakDirection{i}.z;
  angleToAD = angle(peakDirection{i}, vector3d.X) / degree;
  cAxisToADAngleDeg(i) = min(angleToAD, 180 - angleToAD);

  fprintf("OUTPUT=%s MAX_MRD=%.4f C_AXIS_TO_AD_DEG=%.4f\n", ...
    outputFiles(i), maxMRD(i), cAxisToADAngleDeg(i));
end
close all

summary = table(sampleNames, folderNames, inputNames, coldReductionPercent, ...
  totalPixelCount, alphaPixelCount, maxMRD, peakX, peakY, peakZ, ...
  cAxisToADAngleDeg, outputFiles, ...
  'VariableNames', {'sample', 'folder', 'input_file', ...
  'cold_reduction_percent', 'total_pixel_count', 'alpha_ti_pixel_count', ...
  'max_mrd', 'peak_x_ad', 'peak_y_rd', 'peak_z_nd', ...
  'c_axis_to_ad_angle_deg', 'output_file'});

writetable(summary, fullfile(outputDir, "c_axis_pole_figure_summary.csv"));
fprintf("SUMMARY=%s\n", fullfile(outputDir, "c_axis_pole_figure_summary.csv"));
end
