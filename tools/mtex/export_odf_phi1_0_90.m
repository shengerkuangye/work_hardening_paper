function export_odf_phi1_0_90(scanRoot,outputRoot)
%EXPORT_ODF_PHI1_0_90 Evaluate a display window without changing symmetry.
% Each density CSV has rows Phi=0:90 and columns phi1=0:90, in degrees.
arguments
  scanRoot (1,1) string
  outputRoot (1,1) string
end
if ~isfolder(outputRoot), mkdir(outputRoot); end
[summary,fullPeaks,odfs,catalog] = calculate_selected_phi2_odf_data(scanRoot);
writetable(summary,fullfile(outputRoot,"odf_full_space_summary.csv"));
writetable(fullPeaks,fullfile(outputRoot,"odf_full_section_peaks.csv"));
save(fullfile(outputRoot,"odf_models.mat"),"odfs","catalog","summary");
phi2s = [0 30 60];
[x,y] = meshgrid(0:90,0:90);
peaks = fullPeaks;
for i = 1:6
  for j = 1:3
    ori = orientation.byEuler(x(:)*degree,y(:)*degree,phi2s(j)*degree, ...
      "Bunge",odfs{i}.CS,odfs{i}.SS);
    z = reshape(real(eval(odfs{i},ori)),size(x));
    assert(all(isfinite(z),"all"));
    writematrix(z,fullfile(outputRoot,sprintf("density_%s_phi2_%02d.csv", ...
      catalog.sample(i),phi2s(j))));
    [peak,k] = max(z(:));
    r = (i-1)*3+j;
    peaks.section_peak_mrd(r) = peak;
    peaks.phi1_peak_deg(r) = x(k);
    peaks.Phi_peak_deg(r) = y(k);
  end
  fprintf("Exported phi1=0..90 for %s\n",catalog.sample(i));
end
peaks.phi1_min_deg = zeros(18,1);
peaks.phi1_max_deg = repmat(90,18,1);
writetable(peaks,fullfile(outputRoot,"odf_window_peak_positions.csv"));
writematrix(parula(256),fullfile(outputRoot,"parula.csv"));
assert(all(peaks.phi1_peak_deg >= 0 & peaks.phi1_peak_deg <= 90));
assert(all(summary.specimen_symmetry == "1"));
disp("ODF_PHI1_0_90_EXPORT_OK");
end
