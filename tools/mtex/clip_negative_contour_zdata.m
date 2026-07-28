function clippedPointCount = clip_negative_contour_zdata(sectionAxes)
%CLIP_NEGATIVE_CONTOUR_ZDATA Clip nonphysical negative MRD for plotting.

clippedPointCount = 0;
for axisIndex = 1:numel(sectionAxes)
  contourHandles = findall(sectionAxes(axisIndex),"Type","contour");
  assert(numel(contourHandles) == 1, ...
    "Each ODF section axis must contain exactly one contour object.");
  contourValues = contourHandles.ZData;
  assert(all(isfinite(contourValues),"all"), ...
    "ODF contour data must be finite before plotting.");
  negativeMask = contourValues < 0;
  clippedPointCount = clippedPointCount + nnz(negativeMask);
  contourValues(negativeMask) = 0;
  contourHandles.ZData = contourValues;
  assert(all(contourHandles.ZData >= 0,"all"), ...
    "Negative ODF contour data remained after plotting-only clipping.");
end
end
