function [summary, distribution] = generate_c_axis_ad_alignment_metrics( ...
  scanRoot, outputDir)
%GENERATE_C_AXIS_AD_ALIGNMENT_METRICS Analyze alpha-Ti c-axis angles to AD.

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
binEdgesDeg = 0:2:90;

randomMeanAngleDeg = 180 / pi;
randomMedianAngleDeg = 60;
randomAlignmentFactor = 0;
randomFraction15 = 1 - cosd(15);
randomFraction30 = 1 - cosd(30);
randomFraction45 = 1 - cosd(45);

n = numel(sampleNames);
summary = table();
distribution = table();

for i = 1:n
  inputFile = fullfile(scanRoot, folderNames(i), inputNames(i));
  assert(isfile(inputFile), "CTF file not found: %s", inputFile);

  fprintf("IMPORT=%s\n", inputFile);
  ebsd = EBSD.load(inputFile, "convertEuler2SpatialReferenceFrame");
  alpha = ebsd("Ti-Hex");
  assert(~isempty(alpha), "Ti-Hex phase missing in %s", inputFile);

  cAxis = Miller(0, 0, 0, 1, alpha.CS);
  cDirections = alpha.orientations * cAxis;
  thetaDeg = c_axis_angles_to_ad(cDirections);
  [sampleStats, sampleDistribution] = compute_c_axis_ad_statistics( ...
    thetaDeg, length(ebsd), binEdgesDeg);

  summaryRow = struct2table(sampleStats);
  summaryRow = addvars(summaryRow, sampleNames(i), folderNames(i), ...
    inputNames(i), coldReductionPercent(i), 'Before', 1, ...
    'NewVariableNames', {'sample', 'folder', 'input_file', ...
    'cold_reduction_percent'});
  summaryRow = addvars(summaryRow, randomMeanAngleDeg, ...
    randomMedianAngleDeg, randomAlignmentFactor, randomFraction15, ...
    randomFraction30, randomFraction45, ...
    'NewVariableNames', {'random_mean_angle_deg', ...
    'random_median_angle_deg', 'random_alignment_factor', ...
    'random_fraction_within_15deg', 'random_fraction_within_30deg', ...
    'random_fraction_within_45deg'});

  sampleDistribution = addvars(sampleDistribution, ...
    repmat(sampleNames(i), height(sampleDistribution), 1), ...
    repmat(coldReductionPercent(i), height(sampleDistribution), 1), ...
    'Before', 1, 'NewVariableNames', ...
    {'sample', 'cold_reduction_percent'});

  if i == 1
    summary = summaryRow;
    distribution = sampleDistribution;
  else
    summary = [summary; summaryRow]; %#ok<AGROW>
    distribution = [distribution; sampleDistribution]; %#ok<AGROW>
  end

  fprintf("SAMPLE=%s ALPHA_PIXELS=%d MEAN_ANGLE_DEG=%.4f F_AD=%.6f\n", ...
    sampleNames(i), sampleStats.alpha_ti_pixel_count, ...
    sampleStats.mean_angle_deg, sampleStats.alignment_factor);

  clear ebsd alpha cAxis cDirections thetaDeg
end

assert(height(summary) == n, "Expected six summary rows.");
assert(height(distribution) == n * (numel(binEdgesDeg) - 1), ...
  "Unexpected angle-distribution row count.");
assert(all(isfinite(summary{:, 5:end}), "all"), ...
  "Summary contains nonfinite numeric values.");
assert(all(summary.alignment_factor >= -0.5 & ...
  summary.alignment_factor <= 1), "Alignment factor is outside its valid range.");

for i = 1:n
  rows = distribution.sample == sampleNames(i);
  assert(sum(distribution.pixel_count(rows)) == ...
    summary.alpha_ti_pixel_count(i), "Binned pixel count mismatch for %s.", ...
    sampleNames(i));
  assert(abs(sum(distribution.probability(rows)) - 1) < 1e-10, ...
    "Probability sum mismatch for %s.", sampleNames(i));
  assert(abs(sum(distribution.probability_density_per_degree(rows) .* ...
    distribution.bin_width_deg(rows)) - 1) < 1e-10, ...
    "Probability-density integral mismatch for %s.", sampleNames(i));
end

summaryFile = fullfile(outputDir, "c_axis_ad_alignment_summary.csv");
distributionFile = fullfile(outputDir, "c_axis_ad_angle_distribution.csv");
writetable(summary, summaryFile);
writetable(distribution, distributionFile);
fprintf("SUMMARY=%s\n", summaryFile);
fprintf("DISTRIBUTION=%s\n", distributionFile);
end
