function thetaDeg = c_axis_angles_to_ad(cDirections)
%C_AXIS_ANGLES_TO_AD Return antipodal acute angles to specimen x (AD).

assert(isa(cDirections, "vector3d"), ...
  "cDirections must be an MTEX vector3d array.");
thetaDeg = double(angle(cDirections, vector3d.X, "antipodal") / degree);
thetaDeg = thetaDeg(:);
assert(all(isfinite(thetaDeg)) && all(thetaDeg >= 0 & thetaDeg <= 90), ...
  "Computed c-axis angles must be finite and between 0 and 90 degrees.");
end
