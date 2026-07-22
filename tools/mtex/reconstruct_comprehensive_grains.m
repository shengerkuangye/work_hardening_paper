function [grains, ebsdFull, recon] = ...
  reconstruct_comprehensive_grains(ebsdFull, options)
%RECONSTRUCT_COMPREHENSIVE_GRAINS Reconstruct without grid or grain smoothing.

arguments
  ebsdFull EBSD
  options (1,1) struct
end

assert(isfield(options, "detection_threshold_deg"));
assert(isfield(options, "min_grain_pixels"));
detectionThresholdDeg = options.detection_threshold_deg;
minGrainPixels = options.min_grain_pixels;
assert(isscalar(detectionThresholdDeg) && isfinite(detectionThresholdDeg) && ...
  detectionThresholdDeg > 0);
assert(isscalar(minGrainPixels) && isfinite(minGrainPixels) && ...
  minGrainPixels >= 1 && minGrainPixels == fix(minGrainPixels));

sourceIds = double(ebsdFull.id(:));
[grains, grainId] = calcGrains(ebsdFull, 'unitCell', 'threshold', ...
  detectionThresholdDeg * degree);
ebsdFull.grainId = grainId;
assert(isequal(double(ebsdFull.id(:)), sourceIds), ...
  "Grain reconstruction changed persistent EBSD IDs.");

xValues = unique(double(ebsdFull.x(:)));
yValues = unique(double(ebsdFull.y(:)));
outerEndpointIds = double(grains.boundary.ebsdId);
outerEndpointIds = outerEndpointIds(all(outerEndpointIds > 0 & ...
  isfinite(outerEndpointIds), 2), :);
innerEndpointIds = double(grains.innerBoundary.ebsdId);
boundaryEndpointIds = [outerEndpointIds; innerEndpointIds];
nativeGridAudit = audit_native_grid_pairs(boundaryEndpointIds, ...
  sourceIds, [double(ebsdFull.x(:)), double(ebsdFull.y(:))], ...
  median(diff(xValues)), median(diff(yValues)), ...
  "unsmoothed reconstructed grain boundary faces");

summaryMask = grains.numPixel >= minGrainPixels;
recon = struct();
recon.detection_threshold_deg = detectionThresholdDeg;
recon.min_grain_pixels = minGrainPixels;
recon.all_grain_count = length(grains);
recon.summary_grain_count = nnz(summaryMask);
recon.small_grain_count = nnz(~summaryMask);
recon.summary_grain_ids = double(grains.id(summaryMask));
recon.mapped_pixel_count = length(ebsdFull);
recon.indexed_pixel_count = nnz(ebsdFull.isIndexed);
recon.native_grid_audit = nativeGridAudit;
end
