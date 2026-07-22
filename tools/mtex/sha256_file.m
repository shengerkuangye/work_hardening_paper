function digest = sha256_file(filePath)
%SHA256_FILE Return the lowercase SHA-256 digest of a file.

arguments
  filePath (1,1) string
end

assert(isfile(filePath), "File not found: %s", filePath);
fileId = fopen(char(filePath), "rb");
assert(fileId >= 0, "Could not open file for hashing: %s", filePath);
cleanupFile = onCleanup(@() fclose(fileId));
messageDigest = java.security.MessageDigest.getInstance("SHA-256");
chunkBytes = 1024 * 1024;
while true
  bytes = fread(fileId, [chunkBytes, 1], "*uint8");
  if isempty(bytes)
    break
  end
  messageDigest.update(typecast(bytes(:), "int8"));
end
hashBytes = typecast(messageDigest.digest(), "uint8");
digest = string(lower(reshape(dec2hex(hashBytes, 2).', 1, [])));
end
