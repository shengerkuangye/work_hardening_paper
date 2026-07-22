function [inventory, pairs] = generate_comprehensive_ebsd_audit( ...
  scanRoot, outputDir)
%GENERATE_COMPREHENSIVE_EBSD_AUDIT Export scan and raw/denoised audits.

arguments
  scanRoot (1,1) string
  outputDir (1,1) string
end

catalog = comprehensive_ebsd_catalog(scanRoot);
contract = comprehensive_ebsd_output_contract();
auditDir = fullfile(outputDir, contract.directories(1));
if ~isfolder(auditDir)
  mkdir(auditDir);
end

inventoryRows = cell(12, 1);
pairRows = cell(6, 1);
figureHandle = figure("Visible", "off", "Color", "w", ...
  "Position", [20 20 2400 2200]);
cleanupFigure = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 6, 6, "Padding", "compact", ...
  "TileSpacing", "compact");
mapAxes = gobjects(6, 6);
metricMinimum = nan(6, 6);
metricMaximum = nan(6, 6);

for sampleIndex = 1:6
  rawIndex = 2 * sampleIndex - 1;
  denoisedIndex = rawIndex + 1;
  [rawEbsd, rawMeta] = load_comprehensive_ebsd_scan( ...
    catalog(rawIndex, :));
  [denoisedEbsd, denoisedMeta] = load_comprehensive_ebsd_scan( ...
    catalog(denoisedIndex, :));
  inventoryRows{rawIndex} = inventory_row(catalog(rawIndex, :), rawMeta);
  inventoryRows{denoisedIndex} = inventory_row( ...
    catalog(denoisedIndex, :), denoisedMeta);
  [pairRows{sampleIndex}, pixelTable] = compare_raw_denoised_ebsd( ...
    rawEbsd, denoisedEbsd, catalog(rawIndex, :));
  [mapAxes(sampleIndex, :), metricMinimum(sampleIndex, :), ...
    metricMaximum(sampleIndex, :)] = plot_pair_change_maps( ...
    layout, pixelTable, catalog.sample(rawIndex), sampleIndex);
  clear rawEbsd denoisedEbsd pixelTable
end

inventory = vertcat(inventoryRows{:});
pairs = vertcat(pairRows{:});
inventory = inventory(:, cellstr(contract.summaryColumns.scan_inventory));
pairs = pairs(:, cellstr(contract.summaryColumns.raw_denoised_pair_audit));
writetable(inventory, fullfile(auditDir, "scan_inventory.csv"));
writetable(pairs, fullfile(auditDir, "raw_denoised_pair_audit.csv"));
colorLimits = comprehensive_ebsd_change_map_limits( ...
  metricMinimum, metricMaximum);
for metricIndex = 1:6
  for sampleIndex = 1:6
    clim(mapAxes(sampleIndex, metricIndex), colorLimits(metricIndex, :));
  end
end
title(layout, "Raw/denoised EBSD pointwise audit (raw is primary); " + ...
  "differences = denoised - raw");
exportgraphics(figureHandle, ...
  fullfile(auditDir, "raw_denoised_change_maps.png"), ...
  "Resolution", 180);
end

function row = inventory_row(catalogRow, meta)
row = table(string(catalogRow.sample), catalogRow.diameter_mm, ...
  catalogRow.cold_reduction_percent, string(catalogRow.variant), ...
  string(catalogRow.folder), string(catalogRow.input_file), ...
  sha256_file(meta.input_path), meta.x_cells, meta.y_cells, ...
  meta.x_step_um, meta.y_step_um, meta.total_pixels, ...
  meta.indexed_pixels, meta.indexed_fraction, meta.ti_hex_pixels, ...
  meta.ti_cubic_pixels, meta.unindexed_pixels, meta.has_bands, ...
  meta.has_error, meta.has_mad, meta.has_bc, meta.has_bs, ...
  'VariableNames', cellstr(["sample","diameter_mm","cold_reduction_percent", ...
  "variant","folder","input_file","sha256","x_cells","y_cells", ...
  "x_step_um","y_step_um","total_pixels","indexed_pixels", ...
  "indexed_fraction","ti_hex_pixels","ti_cubic_pixels", ...
  "unindexed_pixels","has_bands","has_error","has_mad","has_bc", ...
  "has_bs"]));
end

function [axesRow, minimumRow, maximumRow] = ...
  plot_pair_change_maps(layout, pixelTable, sampleName, rowIndex)
xValues = unique(pixelTable.x_um);
yValues = unique(pixelTable.y_um);
mapValues = { ...
  double(pixelTable.phase_changed), ...
  double(pixelTable.indexing_changed), ...
  pixelTable.orientation_change_deg, ...
  pixelTable.mad_difference, pixelTable.bc_difference, ...
  pixelTable.bs_difference};
metricTitles = ["Phase changed","Indexing changed", ...
  "Orientation change angle (deg)","MAD (denoised - raw)", ...
  "BC (denoised - raw)","BS (denoised - raw)"];
axesRow = gobjects(1, 6);
minimumRow = nan(1, 6);
maximumRow = nan(1, 6);
for columnIndex = 1:6
  axesHandle = nexttile(layout, (rowIndex - 1) * 6 + columnIndex);
  axesRow(columnIndex) = axesHandle;
  imageData = reshape(mapValues{columnIndex}, ...
    [numel(xValues), numel(yValues)]).';
  finiteValues = imageData(isfinite(imageData));
  if ~isempty(finiteValues)
    minimumRow(columnIndex) = min(finiteValues);
    maximumRow(columnIndex) = max(finiteValues);
  end
  imagesc(axesHandle, xValues, yValues, imageData);
  axis(axesHandle, "image");
  set(axesHandle, "YDir", "normal", "FontSize", 7);
  if rowIndex == 1
    title(axesHandle, metricTitles(columnIndex), "FontSize", 9);
  end
  if columnIndex == 1
    ylabel(axesHandle, sampleName + newline + "TD/RD (um)");
  end
  if rowIndex == 6
    xlabel(axesHandle, "AD (um)");
  end
  colorbar(axesHandle);
end
colormap(layout.Parent, turbo(256));
end
