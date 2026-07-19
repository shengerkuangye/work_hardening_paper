function [index1, index2] = four_neighbor_pairs(phaseMask)
%FOUR_NEIGHBOR_PAIRS Return unique horizontal and vertical valid pairs.
assert(islogical(phaseMask) && ismatrix(phaseMask) && ...
  all(size(phaseMask) >= [2,2]), ...
  "phaseMask must be a logical matrix of at least 2-by-2.");
gridIndex = reshape(1:numel(phaseMask), size(phaseMask));

horizontalMask = phaseMask(:,1:end-1) & phaseMask(:,2:end);
horizontal1 = gridIndex(:,1:end-1);
horizontal2 = gridIndex(:,2:end);
verticalMask = phaseMask(1:end-1,:) & phaseMask(2:end,:);
vertical1 = gridIndex(1:end-1,:);
vertical2 = gridIndex(2:end,:);

horizontalIndex1 = horizontal1(horizontalMask);
horizontalIndex2 = horizontal2(horizontalMask);
verticalIndex1 = vertical1(verticalMask);
verticalIndex2 = vertical2(verticalMask);
index1 = [horizontalIndex1(:); verticalIndex1(:)];
index2 = [horizontalIndex2(:); verticalIndex2(:)];
end
