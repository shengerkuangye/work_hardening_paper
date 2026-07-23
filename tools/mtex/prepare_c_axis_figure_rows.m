function prepared = prepare_c_axis_figure_rows(textureSummary)
required = ["sample","diameter_mm","cold_reduction_percent", ...
  "variant","weighting","c_axis_ad_mean_deg", ...
  "c_axis_ad_p10_deg","c_axis_ad_p90_deg"];
assert(istable(textureSummary));
assert(all(ismember(required, ...
  string(textureSummary.Properties.VariableNames))));
textureSummary.sample = string(textureSummary.sample);
textureSummary.variant = string(textureSummary.variant);
textureSummary.weighting = string(textureSummary.weighting);
rows = textureSummary.weighting == "pixel_weighted";
selected = textureSummary(rows, cellstr(required));
assert(height(selected) == 12);
key = selected.sample + "|" + selected.variant;
assert(numel(unique(key)) == 12);
expectedReduction = [0;14.31;26.04;36;43.75;48.98];
expectedSample = ["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"];
expectedDiameter = [7;6.48;6.02;5.6;5.25;5];
prepared = struct();
for variantName = ["raw","denoised"]
  variantRows = selected(selected.variant == variantName,:);
  variantRows = sortrows(variantRows,"cold_reduction_percent");
  assert(height(variantRows) == 6);
  assert(max(abs(variantRows.cold_reduction_percent - ...
    expectedReduction)) < 1e-10);
  assert(isequal(variantRows.sample,expectedSample));
  assert(max(abs(variantRows.diameter_mm - expectedDiameter)) < 1e-10);
  assert(all(isfinite(variantRows{:,6:8}),"all"));
  assert(all(variantRows.c_axis_ad_p10_deg <= ...
    variantRows.c_axis_ad_mean_deg));
  assert(all(variantRows.c_axis_ad_p10_deg >= 0));
  assert(all(variantRows.c_axis_ad_mean_deg <= ...
    variantRows.c_axis_ad_p90_deg));
  assert(all(variantRows.c_axis_ad_p90_deg <= 90));
  prepared.(variantName) = variantRows;
end
prepared.cold_reduction_percent = expectedReduction;
lower = min([prepared.raw.c_axis_ad_p10_deg; ...
  prepared.denoised.c_axis_ad_p10_deg]);
upper = max([prepared.raw.c_axis_ad_p90_deg; ...
  prepared.denoised.c_axis_ad_p90_deg]);
padding = max(1,0.05 * (upper - lower));
prepared.common_y_limits_deg = ...
  [floor(lower-padding),ceil(upper+padding)];
end
