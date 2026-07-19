function [keepMask, componentLength] = component_length_mask( ...
  componentId, segLength, minLength)
%COMPONENT_LENGTH_MASK Retain segments belonging to sufficiently long components.
componentId = double(componentId(:));
segLength = double(segLength(:));
assert(numel(componentId) == numel(segLength));
if isempty(componentId)
  keepMask = false(0,1);
  componentLength = zeros(0,1);
  return
end
assert(all(componentId >= 1 & componentId == fix(componentId)));
assert(all(isfinite(segLength) & segLength > 0));
assert(isscalar(minLength) && isfinite(minLength) && minLength >= 0);
componentLength = accumarray(componentId, segLength, [], @sum);
keepMask = componentLength(componentId) >= minLength;
keepMask = logical(keepMask(:));
end
