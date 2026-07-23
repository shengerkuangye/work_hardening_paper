function test_c_axis_paper_figures(scanRoot, outputRoot)
arguments
  scanRoot (1,1) string = ""
  outputRoot (1,1) string = ""
end

test_prepare_rows();
test_expected_mapping();
test_metadata_builder();
test_input_stats_hash();
if scanRoot ~= ""
  assert(outputRoot ~= "");
  test_formal_generation(scanRoot, outputRoot);
end
fprintf("test_c_axis_paper_figures passed\n");
end

function test_prepare_rows()
sample = repelem(["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"], 4);
diameter_mm = repelem([7;6.48;6.02;5.6;5.25;5], 4);
cold_reduction_percent = repelem([0;14.31;26.04;36;43.75;48.98], 4);
variant = repmat(["raw";"raw";"denoised";"denoised"], 6, 1);
weighting = repmat(["pixel_weighted";"area_weighted_grain_mean"; ...
  "pixel_weighted";"area_weighted_grain_mean"], 6, 1);
c_axis_ad_mean_deg = repelem((76:81)',4) + ...
  repmat([0;20;0.1;20.1],6,1);
c_axis_ad_p10_deg = c_axis_ad_mean_deg - 10;
c_axis_ad_p90_deg = c_axis_ad_mean_deg + 8;
summary = table(sample,diameter_mm,cold_reduction_percent,variant, ...
  weighting,c_axis_ad_mean_deg,c_axis_ad_p10_deg,c_axis_ad_p90_deg);
prepared = prepare_c_axis_figure_rows(summary);
assert(height(prepared.raw) == 6 && height(prepared.denoised) == 6);
assert(all(prepared.raw.weighting == "pixel_weighted"));
assert(isequal(prepared.cold_reduction_percent, ...
  [0;14.31;26.04;36;43.75;48.98]));
assert(all(abs(prepared.denoised.c_axis_ad_mean_deg - ...
  prepared.raw.c_axis_ad_mean_deg - 0.1) < 1e-12));
assert(prepared.common_y_limits_deg(1) <= ...
  min([prepared.raw.c_axis_ad_p10_deg; ...
  prepared.denoised.c_axis_ad_p10_deg]));
assert(prepared.common_y_limits_deg(2) >= ...
  max([prepared.raw.c_axis_ad_p90_deg; ...
  prepared.denoised.c_axis_ad_p90_deg]));

incorrectSample = summary;
incorrectSample.sample(1) = "incorrect";
assert_rejected(@() prepare_c_axis_figure_rows(incorrectSample), ...
  "Incorrect sample mapping must be rejected.");

incorrectDiameter = summary;
incorrectDiameter.diameter_mm(1) = 6.99;
assert_rejected(@() prepare_c_axis_figure_rows(incorrectDiameter), ...
  "Incorrect diameter mapping must be rejected.");

invalidPercentileOrder = summary;
invalidPercentileOrder.c_axis_ad_p10_deg(1) = ...
  invalidPercentileOrder.c_axis_ad_mean_deg(1) + 1;
assert_rejected(@() prepare_c_axis_figure_rows(invalidPercentileOrder), ...
  "P10 greater than the mean must be rejected.");

negativeP10 = summary;
negativeP10.c_axis_ad_p10_deg(1) = -0.01;
assert_rejected(@() prepare_c_axis_figure_rows(negativeP10), ...
  "P10 below zero degrees must be rejected.");

aboveAcuteAngleRange = summary;
aboveAcuteAngleRange.c_axis_ad_p90_deg(1) = 90.01;
assert_rejected(@() prepare_c_axis_figure_rows(aboveAcuteAngleRange), ...
  "P90 above 90 degrees must be rejected.");
end

function test_formal_generation(scanRoot, outputRoot)
catalog = comprehensive_ebsd_catalog(scanRoot);
assert(height(catalog) == 12);
[beforeBytes,beforeTimes,beforeHashes] = ...
  input_stats(catalog.input_path);
metadata = generate_c_axis_paper_figures(scanRoot,outputRoot);
expectedMapping = expected_metadata_mapping();
expectedOutput = expectedMapping.output_name;
assert(isequal(metadata(:,expectedMapping.Properties.VariableNames), ...
  expectedMapping));
assert(height(metadata) == 4);
assert(numel(unique(metadata.pole_color_max_mrd(1:2))) == 1);
assert(numel(unique(metadata.y_min_deg(3:4))) == 1);
assert(numel(unique(metadata.y_max_deg(3:4))) == 1);
textureDir = fullfile(outputRoot,"05_texture");
summary = readtable(fullfile(textureDir,"texture_summary.csv"), ...
  "TextType","string","VariableNamingRule","preserve");
prepared = prepare_c_axis_figure_rows(summary);
expectedTrendRows = {prepared.raw;prepared.denoised};
for trendIndex = 1:2
  metadataIndex = trendIndex + 2;
  expectedRows = expectedTrendRows{trendIndex};
  assert(isequal(metadata.plotted_x{metadataIndex}, ...
    expectedRows.cold_reduction_percent));
  assert(isequal(metadata.plotted_mean_deg{metadataIndex}, ...
    expectedRows.c_axis_ad_mean_deg));
  assert(isequal(metadata.plotted_p10_deg{metadataIndex}, ...
    expectedRows.c_axis_ad_p10_deg));
  assert(isequal(metadata.plotted_p90_deg{metadataIndex}, ...
    expectedRows.c_axis_ad_p90_deg));
end
for fileName = expectedOutput'
  info = dir(fullfile(textureDir,fileName));
  assert(isscalar(info) && info.bytes > 0);
end
[afterBytes,afterTimes,afterHashes] = input_stats(catalog.input_path);
assert(isequal(beforeBytes,afterBytes));
assert(isequal(beforeTimes,afterTimes));
assert(isequal(beforeHashes,afterHashes));
end

function test_expected_mapping()
expectedMapping = expected_metadata_mapping();
assert(height(expectedMapping) == 4);
assert(isequal(expectedMapping.Properties.VariableNames, ...
  {'output_name','variant','figure_kind'}));
end

function test_metadata_builder()
rawPlot = struct();
rawPlot.x = [0;14.31];
rawPlot.mean_deg = [79.1;79.2];
rawPlot.p10_deg = [67.5;69.7];
rawPlot.p90_deg = [88.3;88.2];
denoisedPlot = struct();
denoisedPlot.x = [0;14.31];
denoisedPlot.mean_deg = [79.15;79.25];
denoisedPlot.p10_deg = [67.55;69.75];
denoisedPlot.p90_deg = [88.35;88.25];
metadata = build_c_axis_figure_metadata(4.2,[66 90], ...
  rawPlot,denoisedPlot);
expectedMapping = expected_metadata_mapping();
assert(isequal(metadata(:,expectedMapping.Properties.VariableNames), ...
  expectedMapping));
assert_trend_metadata_matches(metadata,rawPlot,denoisedPlot);

swapped = build_c_axis_figure_metadata(4.2,[66 90], ...
  denoisedPlot,rawPlot);
assert_rejected(@() assert_trend_metadata_matches( ...
  swapped,rawPlot,denoisedPlot), ...
  "Swapped raw and denoised plot data must be rejected.");
end

function assert_trend_metadata_matches(metadata,rawPlot,denoisedPlot)
expectedPlot = {rawPlot;denoisedPlot};
for trendIndex = 1:2
  metadataIndex = trendIndex + 2;
  plotData = expectedPlot{trendIndex};
  assert(isequal(metadata.plotted_x{metadataIndex},plotData.x));
  assert(isequal(metadata.plotted_mean_deg{metadataIndex}, ...
    plotData.mean_deg));
  assert(isequal(metadata.plotted_p10_deg{metadataIndex}, ...
    plotData.p10_deg));
  assert(isequal(metadata.plotted_p90_deg{metadataIndex}, ...
    plotData.p90_deg));
end
end

function expectedMapping = expected_metadata_mapping()
output_name = ["c_axis_pole_figures_raw.png"; ...
  "c_axis_pole_figures_denoised.png"; ...
  "c_axis_mean_orientation_raw.png"; ...
  "c_axis_mean_orientation_denoised.png"];
variant = ["raw";"denoised";"raw";"denoised"];
figure_kind = ["pole_montage";"pole_montage"; ...
  "mean_orientation";"mean_orientation"];
expectedMapping = table(output_name,variant,figure_kind);
end

function test_input_stats_hash()
temporaryPath = string(tempname) + ".txt";
cleanupFile = onCleanup(@() delete(temporaryPath));
fileId = fopen(temporaryPath,"wb");
assert(fileId >= 0);
fwrite(fileId,uint8('abc'),"uint8");
fclose(fileId);
[bytes,~,hashes] = input_stats(temporaryPath);
assert(bytes == 3);
assert(hashes == ...
  "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
clear cleanupFile
end

function [bytes,times,hashes] = input_stats(paths)
bytes = zeros(numel(paths),1);
times = zeros(numel(paths),1);
hashes = strings(numel(paths),1);
for pathIndex = 1:numel(paths)
  info = dir(paths(pathIndex));
  assert(isscalar(info));
  bytes(pathIndex) = info.bytes;
  times(pathIndex) = info.datenum;
  hashes(pathIndex) = sha256_file(paths(pathIndex));
end
end

function assert_rejected(action,message)
didThrow = false;
try
  action();
catch
  didThrow = true;
end
assert(didThrow,message);
end
