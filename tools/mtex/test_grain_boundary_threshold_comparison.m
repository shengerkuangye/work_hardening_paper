function test_grain_boundary_threshold_comparison(scanRoot, outputDir)
arguments
  scanRoot (1,1) string = ""
  outputDir (1,1) string = ""
end

thetaDeg = [2.1;4.9;5;10;14.9;15;30];
segLengthUm = [1;2;3;4;5;6;7];
floorDeg = 2;
indexedAreaUm2 = 100;

stats5 = calculate_boundary_threshold_metrics(thetaDeg, segLengthUm, ...
  floorDeg, 5, indexedAreaUm2);
stats15 = calculate_boundary_threshold_metrics(thetaDeg, segLengthUm, ...
  floorDeg, 15, indexedAreaUm2);

assert(stats5.total_eligible_boundary_length_um == 28);
assert(stats15.total_eligible_boundary_length_um == 28);
assert(stats5.lagb_length_um == 3);
assert(stats15.lagb_length_um == 15);
assert(stats5.hagb_length_um == 25);
assert(stats15.hagb_length_um == 13);
assert(abs(stats5.lagb_length_fraction + stats5.hagb_length_fraction - 1) ...
  < 1e-12);
assert(abs(stats15.lagb_length_fraction + stats15.hagb_length_fraction - 1) ...
  < 1e-12);
assert(stats5.lagb_length_um <= stats15.lagb_length_um);
assert(abs(stats5.lagb_length_density_um_per_um2 - 0.03) < 1e-12);
assert(abs(stats15.lagb_length_density_um_per_um2 - 0.15) < 1e-12);

if scanRoot ~= ""
  assert(outputDir ~= "", "outputDir is required for integration testing.");
  expectedFiles = [
    "grain_boundary_threshold_comparison_summary.csv";
    "grain_boundary_threshold_comparison_distribution.csv"
  ];
  prepare_empty_output_folder(outputDir, expectedFiles);
  generationStart = datetime("now");

  [summary, distribution] = ...
    generate_grain_boundary_threshold_comparison(scanRoot, outputDir);

  expectedSamples = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"];
  expectedThresholds = [5;15];
  assert(height(summary) == 12);
  assert(height(distribution) == 2232);
  assert(isequal(summary.sample, repelem(expectedSamples, 2)));
  assert(isequal(summary.classification_angle_deg, ...
    repmat(expectedThresholds, 6, 1)));
  assert(all(summary.detection_floor_deg == 2));
  assert(all(summary.reconstruction_upper_angle_deg == 15));
  assert(all(summary.min_pixel == 5));
  assert(~any(contains(summary.input_file, "_denoised")));
  assert(all(summary.scan_unit == "um"));
  assert(all(summary.native_x_cell_count == 600));
  assert(all(summary.native_y_cell_count == 600));
  assert(all(summary.nonlocal_ti_ti_endpoint_pair_count == 0));
  assert(all(abs(summary.lagb_length_um + summary.hagb_length_um - ...
    summary.total_eligible_boundary_length_um) < 1e-8));
  assert(all(abs(summary.lagb_length_fraction + ...
    summary.hagb_length_fraction - 1) < 1e-10));

  for sampleName = expectedSamples'
    rows = summary.sample == sampleName;
    row5 = rows & summary.classification_angle_deg == 5;
    row15 = rows & summary.classification_angle_deg == 15;
    assert(nnz(row5) == 1 && nnz(row15) == 1);
    assert(abs(summary.total_eligible_boundary_length_um(row5) - ...
      summary.total_eligible_boundary_length_um(row15)) < 1e-8);
    assert(summary.lagb_length_um(row5) <= ...
      summary.lagb_length_um(row15) + 1e-8);

    distribution5 = distribution(distribution.sample == sampleName & ...
      distribution.classification_angle_deg == 5, :);
    distribution15 = distribution(distribution.sample == sampleName & ...
      distribution.classification_angle_deg == 15, :);
    assert(height(distribution5) == 186 && height(distribution15) == 186);
    assert(isequal(distribution5.bin_center_deg, ...
      distribution15.bin_center_deg));
    assert(all(abs(distribution5.boundary_length_um - ...
      distribution15.boundary_length_um) < 1e-8));
    assert(abs(sum(distribution5.length_fraction) - 1) < 1e-10);
    assert(abs(sum(distribution15.length_fraction) - 1) < 1e-10);
  end

  generatorText = string(fileread(which( ...
    "generate_grain_boundary_threshold_comparison")));
  assert(contains(generatorText, "calcGrains(ebsdFull, 'unitCell'"));
  assert(contains(generatorText, "'threshold', [2 15] * degree"));
  forbiddenCalls = [
    "calcGrains(ebsd('indexed')";
    "calcGrains(ebsd(""indexed"")";
    "calcGrains(ebsdFull('indexed')";
    "calcGrains(ebsdFull(""indexed"")";
    "gridify(";
    "interpolate(";
    "smooth("
  ];
  assert(~any(contains(generatorText, forbiddenCalls)));

  summaryRead = readtable(fullfile(outputDir, ...
    "grain_boundary_threshold_comparison_summary.csv"), ...
    "TextType", "string");
  distributionRead = readtable(fullfile(outputDir, ...
    "grain_boundary_threshold_comparison_distribution.csv"), ...
    "TextType", "string");
  assert_table_round_trip(summaryRead, summary);
  assert_table_round_trip(distributionRead, distribution);
  assert_fresh_exact_artifacts(outputDir, expectedFiles, generationStart);
end

fprintf("GRAIN_BOUNDARY_THRESHOLD_COMPARISON_TESTS_OK\n");
end

function assert_table_round_trip(actual, expected)
assert(isequal(actual.Properties.VariableNames, ...
  expected.Properties.VariableNames));
assert(height(actual) == height(expected));
for variableName = string(expected.Properties.VariableNames)
  actualValues = actual.(variableName);
  expectedValues = expected.(variableName);
  if isnumeric(expectedValues)
    assert(isequal(size(actualValues), size(expectedValues)));
    finiteMask = isfinite(actualValues) & isfinite(expectedValues);
    assert(all(isnan(actualValues(~finiteMask)) == ...
      isnan(expectedValues(~finiteMask))));
    assert(all(abs(actualValues(finiteMask) - ...
      expectedValues(finiteMask)) < 1e-10));
  else
    assert(isequal(string(actualValues), string(expectedValues)));
  end
end
end

function prepare_empty_output_folder(outputDir, expectedFiles)
if ~isfolder(outputDir)
  mkdir(outputDir);
end
for fileName = expectedFiles'
  filePath = fullfile(outputDir, fileName);
  if isfile(filePath)
    delete(filePath);
  end
end
remaining = dir(outputDir);
remaining = remaining(~[remaining.isdir]);
assert(isempty(remaining), ...
  "Integration output folder contains an unrelated file.");
end

function assert_fresh_exact_artifacts(outputDir, expectedFiles, generationStart)
files = dir(outputDir);
files = files(~[files.isdir]);
actualNames = sort(string({files.name}));
expectedNames = sort(expectedFiles).';
assert(isequal(actualNames, expectedNames));
for fileIndex = 1:numel(files)
  assert(files(fileIndex).bytes > 0);
  assert(datetime(files(fileIndex).datenum, "ConvertFrom", "datenum") >= ...
    generationStart - seconds(2));
end
end
