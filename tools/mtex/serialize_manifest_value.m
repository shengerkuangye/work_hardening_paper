function result = serialize_manifest_value(value)
%SERIALIZE_MANIFEST_VALUE Convert one registered value to one CSV cell.

if isstring(value)
  if isscalar(value)
    result = value;
  else
    result = string(jsonencode(value));
  end
elseif ischar(value)
  result = string(value);
elseif isnumeric(value) || islogical(value)
  result = string(mat2str(value));
else
  result = string(jsonencode(value));
end
assert(isscalar(result));
end
