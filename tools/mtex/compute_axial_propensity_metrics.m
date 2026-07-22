function [byGrain, summary, meta] = compute_axial_propensity_metrics( ...
  orientations, grainIds, grainAreas, options)
%COMPUTE_AXIAL_PROPENSITY_METRICS Geometric HCP propensity under AD tension.

arguments
  orientations orientation
  grainIds (:,1) double
  grainAreas (:,1) double
  options (1,1) struct = struct()
end

assert(~isempty(orientations), "orientations must not be empty.");
nGrains = length(orientations);
assert(numel(grainIds) == nGrains && numel(grainAreas) == nGrains, ...
  "Orientation, grain-ID, and area counts must match.");
assert(all(isfinite(grainIds) & grainIds > 0 & ...
  grainIds == fix(grainIds)) && numel(unique(grainIds)) == nGrains, ...
  "grainIds must be unique positive integers.");
assert(all(isfinite(grainAreas) & grainAreas > 0), ...
  "grainAreas must contain finite positive values.");

includeTwins = option_or_default(options, "include_twins", true);
assert(isscalar(includeTwins) && (islogical(includeTwins) || ...
  isnumeric(includeTwins)), "include_twins must be scalar logical.");
includeTwins = logical(includeTwins);
taylorCrss = option_or_default(options, "taylor_crss", []);
if ~isempty(taylorCrss)
  taylorCrss = double(taylorCrss(:)');
  assert(numel(taylorCrss) == 4 && ...
    all(isfinite(taylorCrss) & taylorCrss > 0), ...
    "taylor_crss must be empty or four finite positive values.");
end

cs = orientations.CS;
familyNames = ["basal_a"; "prismatic_a"; "pyramidal_a"; ...
  "pyramidal_ca"];
familySystems = {slipSystem.basal(cs), ...
  slipSystem.prismaticA(cs), slipSystem.pyramidalA(cs), ...
  slipSystem.pyramidalCA(cs)};
isTwinFamily = false(4, 1);
if includeTwins
  familyNames = [familyNames; "extension_twin_t1"; ...
    "contraction_twin_c1"];
  familySystems = [familySystems, {slipSystem.twinT1(cs), ...
    slipSystem.twinC1(cs)}];
  isTwinFamily = [isTwinFamily; true; true];
end

[taylorFactor, taylorStatus, taylorText] = calculate_taylor_sensitivity( ...
  orientations, cs, taylorCrss);
areaWeight = grainAreas / sum(grainAreas);
[phi1, Phi, phi2] = Euler(orientations);
phi1Deg = double(phi1(:) / degree);
PhiDeg = double(Phi(:) / degree);
phi2Deg = double(phi2(:) / degree);

nFamilies = numel(familyNames);
nRows = nGrains * nFamilies;
grain_id = zeros(nRows, 1);
family = strings(nRows, 1);
max_abs_schmid = zeros(nRows, 1);
winning_variant = zeros(nRows, 1);
area_um2 = zeros(nRows, 1);
area_weight = zeros(nRows, 1);
phi1_deg = zeros(nRows, 1);
Phi_deg = zeros(nRows, 1);
phi2_deg = zeros(nRows, 1);
taylor_factor = nan(nRows, 1);
variantCounts = zeros(nFamilies, 1);

for familyIndex = 1:nFamilies
  variants = familySystems{familyIndex}.symmetrise("antipodal");
  variantCounts(familyIndex) = length(variants);
  specimenSystems = orientations * variants;
  schmid = reshape(abs(specimenSystems.SchmidFactor(vector3d.X)), ...
    nGrains, []);
  assert(all(isfinite(schmid), "all"), ...
    "Schmid-factor calculation returned nonfinite values.");
  assert(all(schmid >= -1e-12 & schmid <= 0.5 + 1e-10, "all"), ...
    "Absolute Schmid factors must lie in [0, 0.5].");
  [familyMax, winner] = max(schmid, [], 2);
  rows = (familyIndex - 1) * nGrains + (1:nGrains);
  grain_id(rows) = grainIds;
  family(rows) = familyNames(familyIndex);
  max_abs_schmid(rows) = familyMax;
  winning_variant(rows) = winner;
  area_um2(rows) = grainAreas;
  area_weight(rows) = areaWeight;
  phi1_deg(rows) = phi1Deg;
  Phi_deg(rows) = PhiDeg;
  phi2_deg(rows) = phi2Deg;
  taylor_factor(rows) = taylorFactor;
end

byGrain = table(grain_id, family, max_abs_schmid, winning_variant, ...
  area_um2, area_weight, phi1_deg, Phi_deg, phi2_deg, taylor_factor);

assumption = strings(nFamilies, 1);
grain_count = repmat(nGrains, nFamilies, 1);
area_weighted_mean_max_schmid = zeros(nFamilies, 1);
area_weighted_median_max_schmid = zeros(nFamilies, 1);
area_fraction_schmid_ge_0_4 = zeros(nFamilies, 1);
taylor_factor_area_weighted_mean = nan(nFamilies, 1);
for familyIndex = 1:nFamilies
  rows = byGrain.family == familyNames(familyIndex);
  values = byGrain.max_abs_schmid(rows);
  weights = byGrain.area_weight(rows);
  assumption(familyIndex) = ...
    "geometric absolute Schmid factor; no CRSS activity claim; " + ...
    taylorText;
  if isTwinFamily(familyIndex)
    assumption(familyIndex) = assumption(familyIndex) + ...
      "; absolute-value twin screen does not resolve polarity";
  end
  area_weighted_mean_max_schmid(familyIndex) = sum(weights .* values);
  area_weighted_median_max_schmid(familyIndex) = ...
    weighted_median(values, weights);
  area_fraction_schmid_ge_0_4(familyIndex) = ...
    sum(weights(values >= 0.4));
  if all(isfinite(taylorFactor))
    taylor_factor_area_weighted_mean(familyIndex) = ...
      sum(areaWeight .* taylorFactor);
  end
end
summary = table(familyNames, assumption, grain_count, ...
  area_weighted_mean_max_schmid, ...
  area_weighted_median_max_schmid, ...
  area_fraction_schmid_ge_0_4, ...
  taylor_factor_area_weighted_mean, ...
  'VariableNames', {'family', 'assumption', 'grain_count', ...
  'area_weighted_mean_max_schmid', ...
  'area_weighted_median_max_schmid', ...
  'area_fraction_schmid_ge_0_4', ...
  'taylor_factor_area_weighted_mean'});

meta = struct();
meta.load_axis = "AD";
meta.load_axis_vector = "vector3d.X";
meta.family_names = familyNames;
meta.family_variant_count = variantCounts;
meta.taylor_status = taylorStatus;
meta.taylor_crss = taylorCrss;
end

function value = option_or_default(options, fieldName, defaultValue)
if isfield(options, fieldName)
  value = options.(fieldName);
else
  value = defaultValue;
end
end

function [taylorFactor, status, text] = ...
  calculate_taylor_sensitivity(orientations, cs, crss)
nGrains = length(orientations);
taylorFactor = nan(nGrains, 1);
if isempty(crss)
  status = "not_requested_no_crss";
  text = "Taylor sensitivity not computed because CRSS was not supplied";
  return
end

baseSystems = [slipSystem.basal(cs, crss(1)), ...
  slipSystem.prismaticA(cs, crss(2)), ...
  slipSystem.pyramidalA(cs, crss(3)), ...
  slipSystem.pyramidalCA(cs, crss(4))];
systems = baseSystems.symmetrise;
specimenStrain = strainTensor(diag([1 -0.5 -0.5]));
try
  calculated = double(calcTaylor(inv(orientations) * ...
    specimenStrain, systems)); %#ok<MINV>
  calculated = calculated(:);
  if numel(calculated) ~= nGrains || ...
      any(~isfinite(calculated) | calculated <= 0)
    status = "failed_nonfinite_or_nonpositive_result";
    text = "Taylor sensitivity failed; values exported as NaN";
    return
  end
  taylorFactor = calculated;
  details = "combined slip set; family/CRSS order " + ...
    "basal_a,prismatic_a,pyramidal_a,pyramidal_ca; " + ...
    "eps_AD=diag([1,-0.5,-0.5]); " + ...
    "CRSS are assumed relative weights; " + ...
    "same combined Taylor sensitivity is repeated on every family row";
  if max(crss) - min(crss) < 10 * eps(max(crss))
    status = "ok_equal_crss";
    text = "Taylor sensitivity equal CRSS=[" + ...
      join(string(compose("%.6g", crss)), ",") + "]; " + ...
      details;
  else
    status = "ok_supplied_crss";
    text = "Taylor sensitivity supplied CRSS=[" + ...
      join(string(compose("%.6g", crss)), ",") + "]; " + ...
      details;
  end
catch exception
  status = "failed_" + string(exception.identifier);
  text = "Taylor sensitivity failed (" + ...
    string(exception.identifier) + "); values exported as NaN";
end
end

function result = weighted_median(values, weights)
[sortedValues, order] = sort(values);
sortedWeights = weights(order) / sum(weights);
result = sortedValues(find(cumsum(sortedWeights) >= 0.5, 1, "first"));
end
