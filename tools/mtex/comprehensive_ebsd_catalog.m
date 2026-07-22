function catalog = comprehensive_ebsd_catalog(scanRoot)
%COMPREHENSIVE_EBSD_CATALOG Register raw and denoised EBSD inputs.

arguments
  scanRoot (1,1) string
end

assert(isfolder(scanRoot), "EBSD scan folder not found: %s", scanRoot);

sample = repelem(["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"], 2);
diameter_mm = repelem([7;6.48;6.02;5.6;5.25;5], 2);
cold_reduction_percent = repelem( ...
  [0;14.31;26.04;36.00;43.75;48.98], 2);
variant = repmat(["raw";"denoised"], 6, 1);
folder = repelem(["d7";"d6_48";"d6_02";"d5_6";"d5_25";"d5"], 2);
rawNames = [ ...
  "ebsd_sample_7_map_15.ctf"
  "ebsd_sample_648_map_13.ctf"
  "ebsd_sample_602_map_11.ctf"
  "ebsd_sample_56_map_9.ctf"
  "ebsd_sample_525_map_7.ctf"
  "ebsd_sample_5_map_3.ctf"
];
input_file = strings(12, 1);
input_file(1:2:end) = rawNames;
input_file(2:2:end) = replace(rawNames, ".ctf", "_denoised.ctf");

rootPath = canonical_path(scanRoot);
input_path = strings(12, 1);
for rowIndex = 1:12
  input_path(rowIndex) = canonical_path(fullfile(rootPath, ...
    folder(rowIndex), input_file(rowIndex)));
end
coordinate_x = repmat("AD", 12, 1);
coordinate_y = repmat("TD_RD", 12, 1);
coordinate_z = repmat("ND", 12, 1);
x_cells = repmat(600, 12, 1);
y_cells = repmat(600, 12, 1);
x_step_um = repmat(0.5, 12, 1);
y_step_um = repmat(0.5, 12, 1);

catalog = table(sample, diameter_mm, cold_reduction_percent, variant, ...
  folder, input_file, input_path, coordinate_x, coordinate_y, ...
  coordinate_z, x_cells, y_cells, x_step_um, y_step_um);
assert(all(arrayfun(@isfile, catalog.input_path)), ...
  "One or more registered CTF files do not exist.");
assert(numel(unique(catalog.input_path)) == height(catalog), ...
  "The EBSD catalog contains duplicate input paths.");
end

function path = canonical_path(path)
path = string(char(java.io.File(char(path)).getCanonicalPath()));
end
