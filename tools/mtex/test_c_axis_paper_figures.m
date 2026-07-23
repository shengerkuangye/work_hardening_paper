function test_c_axis_paper_figures(scanRoot, outputRoot)
arguments
  scanRoot (1,1) string = ""
  outputRoot (1,1) string = ""
end

test_prepare_rows();
if scanRoot ~= ""
  assert(outputRoot ~= "");
  test_formal_generation(scanRoot, outputRoot);
end
fprintf("test_c_axis_paper_figures passed\n");
end

function test_prepare_rows()
sample = repelem(["7d";"6.48d";"6.02d";"5.6d";"5.25d";"5d"], 4);
diameter_mm = repelem([7;6.48;6.02;5.6;5.25;5], 4);
cold_reduction_percent = repelem([0;14.31;26.04;36;43.75;48.98], 4);
variant = repmat(["raw";"raw";"denoised";"denoised"], 6, 1);
weighting = repmat(["pixel_weighted";"area_weighted_grain_mean"; ...
  "pixel_weighted";"area_weighted_grain_mean"], 6, 1);
c_axis_ad_mean_deg = repelem((79:84)',4) + ...
  repmat([0;20;0.1;20.1],6,1);
c_axis_ad_p10_deg = c_axis_ad_mean_deg - 10;
c_axis_ad_p90_deg = c_axis_ad_mean_deg + 8;
summary = table(sample,diameter_mm,cold_reduction_percent,variant, ...
  weighting,c_axis_ad_mean_deg,c_axis_ad_p10_deg,c_axis_ad_p90_deg);
prepared = prepare_c_axis_figure_rows(summary);
assert(height(prepared.raw) == 6 && height(prepared.denoised) == 6);
assert(all(prepared.raw.weighting == "pixel_weighted"));
assert(isequal(prepared.cold_reduction_percent, ...
  [0;14.31;26.04;36;43.75;48.98]));
assert(all(abs(prepared.denoised.c_axis_ad_mean_deg - ...
  prepared.raw.c_axis_ad_mean_deg - 0.1) < 1e-12));
assert(prepared.common_y_limits_deg(1) <= ...
  min([prepared.raw.c_axis_ad_p10_deg; ...
  prepared.denoised.c_axis_ad_p10_deg]));
assert(prepared.common_y_limits_deg(2) >= ...
  max([prepared.raw.c_axis_ad_p90_deg; ...
  prepared.denoised.c_axis_ad_p90_deg]));
end
