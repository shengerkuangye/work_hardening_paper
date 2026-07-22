function [filteredEndpointIds, removedRows] = ...
  validate_filter_boundary_endpoint_ids(endpointIds, boundaryKind)
%VALIDATE_FILTER_BOUNDARY_ENDPOINT_IDS Validate native boundary point IDs.

arguments
  endpointIds (:,2) double
  boundaryKind (1,1) string {mustBeMember(boundaryKind, ...
    ["outer","inner"])}
end

invalid = ~isfinite(endpointIds) | endpointIds < 0 | ...
  endpointIds ~= fix(endpointIds);
if any(invalid, "all")
  error("validate_filter_boundary_endpoint_ids:InvalidEndpointIds", ...
    "Boundary endpoint IDs must be nonnegative finite integers; " + ...
    "invalid values cannot be discarded.");
end

zeroCount = sum(endpointIds == 0, 2);
if boundaryKind == "outer" && any(zeroCount == 2)
  error("validate_filter_boundary_endpoint_ids:DegenerateOuterEndpoint", ...
    "An outer scan-perimeter boundary row must contain exactly one " + ...
    "zero endpoint and one positive endpoint.");
end
if boundaryKind == "inner" && any(zeroCount > 0)
  error("validate_filter_boundary_endpoint_ids:UnexpectedZeroEndpoint", ...
    "Only outer scan-perimeter boundary rows may contain zero endpoints.");
end
removedRows = boundaryKind == "outer" & zeroCount == 1;
filteredEndpointIds = endpointIds(~removedRows, :);
end
