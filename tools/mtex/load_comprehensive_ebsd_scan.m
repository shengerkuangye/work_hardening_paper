function [ebsdFull, meta] = load_comprehensive_ebsd_scan(catalogRow)
%LOAD_COMPREHENSIVE_EBSD_SCAN Load and validate one registered native grid.

assert(istable(catalogRow) && height(catalogRow) == 1, ...
  "catalogRow must be a one-row table from comprehensive_ebsd_catalog.");
assert(~isempty(which("EBSD")), ...
  "MTEX is not loaded in the current MATLAB session.");
inputPath = string(catalogRow.input_path);
assert(isfile(inputPath), "CTF file not found: %s", inputPath);

ebsdFull = EBSD.load(inputPath, "convertEuler2SpatialReferenceFrame");
xValues = unique(double(ebsdFull.x(:)));
yValues = unique(double(ebsdFull.y(:)));
xSteps = diff(xValues);
ySteps = diff(yValues);
coordinateTolerance = 1e-10;

validGeometry = string(ebsdFull.scanUnit) == "um" && ...
  numel(xValues) == catalogRow.x_cells && ...
  numel(yValues) == catalogRow.y_cells && ...
  length(ebsdFull) == catalogRow.x_cells * catalogRow.y_cells && ...
  all(abs(xSteps - catalogRow.x_step_um) < coordinateTolerance) && ...
  all(abs(ySteps - catalogRow.y_step_um) < coordinateTolerance) && ...
  size(unique([double(ebsdFull.x(:)), double(ebsdFull.y(:))], ...
  "rows"), 1) == length(ebsdFull);
if ~validGeometry
  error("load_comprehensive_ebsd_scan:GeometryMismatch", ...
    "Registered geometry does not match %s.", inputPath);
end

ids = double(ebsdFull.id(:));
assert(numel(unique(ids)) == length(ebsdFull) && ...
  all(isfinite(ids) & ids > 0 & ids == fix(ids)), ...
  "Persistent EBSD IDs must be unique positive integers in %s.", ...
  inputPath);

qualityFields = string(fieldnames(ebsdFull.prop));
indexedPixels = nnz(ebsdFull.isIndexed);
meta = struct();
meta.input_path = inputPath;
meta.mtex_version = string(getMTEXpref("version"));
meta.coordinate_x = string(catalogRow.coordinate_x);
meta.coordinate_y = string(catalogRow.coordinate_y);
meta.coordinate_z = string(catalogRow.coordinate_z);
meta.x_cells = numel(xValues);
meta.y_cells = numel(yValues);
meta.x_step_um = median(xSteps);
meta.y_step_um = median(ySteps);
meta.total_pixels = length(ebsdFull);
meta.indexed_pixels = indexedPixels;
meta.indexed_fraction = indexedPixels / meta.total_pixels;
meta.ti_hex_pixels = length(ebsdFull("Ti-Hex"));
meta.ti_cubic_pixels = length(ebsdFull("Titanium cubic"));
meta.unindexed_pixels = length(ebsdFull("notIndexed"));
meta.phase_counts = struct("ti_hex", meta.ti_hex_pixels, ...
  "ti_cubic", meta.ti_cubic_pixels, ...
  "unindexed", meta.unindexed_pixels);
meta.has_bands = ismember("bands", qualityFields);
meta.has_error = ismember("error", qualityFields);
meta.has_mad = ismember("mad", qualityFields);
meta.has_bc = ismember("bc", qualityFields);
meta.has_bs = ismember("bs", qualityFields);
end
