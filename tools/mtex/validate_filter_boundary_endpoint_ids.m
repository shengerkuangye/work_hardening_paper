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

zeroRows = any(endpointIds == 0, 2);
if boundaryKind == "inner" && any(zeroRows)
  error("validate_filter_boundary_endpoint_ids:UnexpectedZeroEndpoint", ...
    "Only outer scan-perimeter boundary rows may contain zero endpoints.");
end
removedRows = boundaryKind == "outer" & zeroRows;
filteredEndpointIds = endpointIds(~removedRows, :);
end
