# Comprehensive EBSD Analysis Implementation Plan

> **Required sub-skill:** Use `superpowers:subagent-driven-development` to execute this plan. After Task 2 establishes the shared interfaces, use `superpowers:dispatching-parallel-agents` for Tasks 3-5.

**Goal:** Build a reproducible MTEX 6.1.1 workflow that extracts the maximum defensible information from all six original CTF scans, uses the six denoised CTF scans only as paired sensitivity evidence, and integrates EBSD trends with axial tensile behavior without overstating causality.

**Architecture:** A tested shared foundation owns sample metadata, coordinate conventions, import, grain reconstruction, and audit fields. Three independent analysis packages then produce maps/morphology/boundaries, intragranular deformation/texture, and axial slip/twinning propensity. A final runner merges their machine-readable outputs with tensile summaries and creates manuscript-oriented comparison figures and a provenance manifest.

**Tech stack:** MATLAB R2025a, MTEX 6.1.1, MATLAB `matlab.unittest`-style assertion functions, CSV/PNG/SVG outputs, PowerShell SHA-256 verification.

## Global constraints

- Never modify files below `data/` or `references/`.
- Treat the map horizontal direction as the bar axial direction: specimen `x = AD`; map vertical direction `y = TD/RD`; specimen `z = ND`.
- Import with `EBSD.load(file,'convertEuler2SpatialReferenceFrame')`; preserve the native 600 x 600 grid and 0.5 micrometre step.
- Original CTF values are the primary quantitative results. Denoised values are paired robustness checks and may be used for cleaner illustrative maps only when labeled.
- Reconstruct grains without geometric smoothing. Use a 2 degree detection floor, 15 degree LAGB/HAGB classification, and `minPixel = 5` as the registered primary analysis; report 0.5, 1, 2, and 5 degree sensitivity where applicable.
- Weight boundary statistics by physical segment length; report count-weighted values only as sensitivity results.
- Do not call KAM, GOS, GROD, texture intensity, Schmid factor, or Taylor factor a direct dislocation-density measurement.
- Schmid/Taylor results describe orientation-based propensity under a subsequent uniaxial tensile load parallel to AD; they do not reconstruct the multiaxial rotary-swaging stress path.
- Every generated CSV includes `sample`, `diameter_mm`, `cold_reduction_percent`, `variant`, and the analysis parameters relevant to that row.
- Every figure labels raw/denoised provenance and uses identical limits within a comparison panel.

## Execution graph

```text
Task 1 tests and output contract
          |
Task 2 shared foundation + audit
     /          |           \
Task 3       Task 4       Task 5       (parallel, disjoint files)
     \          |           /
          Task 6 integration
                 |
          Task 7 full verification
```

### Task 1: Lock the output contract with failing tests

**Files:**

- Create: `tools/mtex/test_comprehensive_ebsd_contract.m`
- Create: `tools/mtex/comprehensive_ebsd_output_contract.m`

**Step 1: Write the contract test first**

The test must assert that the contract returns the following eight output directories and required artifacts:

- `00_audit`: `scan_inventory.csv`, `raw_denoised_pair_audit.csv`, `raw_denoised_change_maps.png`
- `01_standard_maps`: one `sample_variant_maps.png` per scan variant
- `02_grain_morphology`: `grain_morphology_by_grain.csv`, `grain_morphology_summary.csv`, `grain_morphology_trends.png`
- `03_boundaries`: `boundary_segments.csv`, `boundary_summary.csv`, `boundary_angle_distributions.csv`, `boundary_trends.png`
- `04_intragranular`: `intragranular_by_pixel.csv`, `intragranular_by_grain.csv`, `intragranular_summary.csv`, `intragranular_trends.png`
- `05_texture`: `texture_summary.csv`, `c_axis_orientation_distribution.csv`, `pole_figures.png`, `inverse_pole_figures.png`, `odf_sections.png`, `texture_trends.png`
- `06_axial_propensity`: `axial_propensity_by_grain.csv`, `axial_propensity_summary.csv`, `axial_propensity_trends.png`
- `07_raw_denoised_comparison`: `raw_denoised_metric_comparison.csv`, `raw_denoised_metric_comparison.png`
- `08_tensile_integration`: `ebsd_tensile_merged.csv`, `ebsd_tensile_rank_correlations.csv`, `ebsd_tensile_comparison.png`
- root: `README.md`, `analysis_manifest.csv`

The contract must also expose the registered parameters and the exact column sets for all summary tables.

**Step 2: Run the test and observe failure**

Run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath('tools/mtex'); test_comprehensive_ebsd_contract"
```

Expected: failure because `comprehensive_ebsd_output_contract` does not yet exist.

**Step 3: Implement the smallest contract function**

Return a scalar struct with `directories`, `artifacts`, `summaryColumns`, and `parameters`. Avoid embedding analysis logic.

**Step 4: Re-run the test**

Expected: `test_comprehensive_ebsd_contract passed` and exit code 0.

**Step 5: Commit**

```powershell
git add tools/mtex/test_comprehensive_ebsd_contract.m tools/mtex/comprehensive_ebsd_output_contract.m
git commit -m "test: define comprehensive EBSD output contract"
```

### Task 2: Build the shared catalog, loader, reconstruction, and paired audit

**Files:**

- Create: `tools/mtex/comprehensive_ebsd_catalog.m`
- Create: `tools/mtex/load_comprehensive_ebsd_scan.m`
- Create: `tools/mtex/reconstruct_comprehensive_grains.m`
- Create: `tools/mtex/compare_raw_denoised_ebsd.m`
- Create: `tools/mtex/generate_comprehensive_ebsd_audit.m`
- Create: `tools/mtex/test_comprehensive_ebsd_foundation.m`
- Modify: `tools/mtex/test_comprehensive_ebsd_contract.m`

**Interfaces:**

```matlab
catalog = comprehensive_ebsd_catalog(scanRoot)
[ebsdFull, meta] = load_comprehensive_ebsd_scan(catalogRow)
[grains, ebsdFull, recon] = reconstruct_comprehensive_grains(ebsdFull, options)
[pairRow, pixelTable] = compare_raw_denoised_ebsd(rawEbsd, denoisedEbsd, sampleMeta)
[inventory, pairs] = generate_comprehensive_ebsd_audit(scanRoot, outputDir)
```

`catalog` has 12 rows in deformation order: diameters 7, 6.48, 6.02, 5.6, 5.25, and 5 mm, each followed by raw and denoised variants. The associated cold area reductions are 0, 14.31, 26.04, 36.00, 43.75, and 48.98 percent.

**Step 1: Write synthetic and integration tests**

Test catalog ordering, file existence, absence of duplicated paths, correct AD/TD/ND labels, the 600 x 600 grid, 0.5 micrometre step, and identical raw/denoised coordinates. Assert that paired comparison conserves pixel count and reports phase changes, indexing changes, Euler/orientation change angle, MAD, BC, BS, Bands, and Error differences.

**Step 2: Run the foundation test and observe failure**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath('C:/Users/22069/Documents/MATLAB/mtex-6.1.1'); startup_mtex; addpath('tools/mtex'); test_comprehensive_ebsd_foundation('data/ebsd_kpl_250221_7_df/scans')"
```

Expected: missing-function failure.

**Step 3: Implement catalog and loader**

The loader records the absolute input path, MTEX version, grid dimensions, step, phase counts, indexed fraction, and quality-field availability. It rejects any scan whose coordinates or dimensions violate the registered geometry.

**Step 4: Implement unsmoothed grain reconstruction**

Use `calcGrains` on the full mapped scan with the requested detection threshold; retain persistent EBSD IDs. Filter summary grains with `numPixel >= 5`, but do not erase small grains from audit totals. Reuse `audit_native_grid_pairs.m` for neighbor-network integrity.

**Step 5: Implement paired audit and figures**

For each pair, export a single row of summary differences and pixel-level spatial maps for changed phase/indexing/orientation plus MAD/BC/BS change. Do not overwrite either CTF.

**Step 6: Re-run tests**

Expected: all 12 scans load, all six raw/denoised pairs align exactly by native coordinate, and the audit artifacts match the contract.

**Step 7: Commit**

```powershell
git add tools/mtex/comprehensive_ebsd_catalog.m tools/mtex/load_comprehensive_ebsd_scan.m tools/mtex/reconstruct_comprehensive_grains.m tools/mtex/compare_raw_denoised_ebsd.m tools/mtex/generate_comprehensive_ebsd_audit.m tools/mtex/test_comprehensive_ebsd_foundation.m tools/mtex/test_comprehensive_ebsd_contract.m
git commit -m "feat: add EBSD foundation and raw-denoised audit"
```

### Task 3: Generate standard maps, grain morphology, and boundary statistics

**Files:**

- Create: `tools/mtex/compute_grain_morphology_metrics.m`
- Create: `tools/mtex/compute_boundary_network_metrics.m`
- Create: `tools/mtex/generate_comprehensive_maps_morphology_boundaries.m`
- Create: `tools/mtex/test_comprehensive_maps_morphology_boundaries.m`
- Reuse: `tools/mtex/partition_ti_hex_boundary_segments.m`
- Reuse: `tools/mtex/calculate_boundary_threshold_metrics.m`
- Reuse: `tools/mtex/audit_native_grid_pairs.m`

**Step 1: Write synthetic geometry and boundary tests**

Use synthetic grains to test equivalent diameter, area, perimeter, `fitEllipse` long/short axes, aspect ratio, `longAxis` angle folded to 0-90 degrees from AD, `caliper` maximum/minimum Feret diameter, and `shapeFactor`. Use synthetic boundary segments to test length conservation and bins 2-5, 5-15, and at least 15 degrees.

**Step 2: Run and observe failure**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath('C:/Users/22069/Documents/MATLAB/mtex-6.1.1'); startup_mtex; addpath('tools/mtex'); test_comprehensive_maps_morphology_boundaries"
```

**Step 3: Implement morphology tables**

Create per-grain and per-scan summaries with grain-count and area-weighted medians/quantiles. Include boundary-touching flags and provide a summary sensitivity excluding boundary-touching grains.

**Step 4: Implement standard map panels**

For each variant export IPF-X/AD, IPF-Y/TD, IPF-Z/ND, BC, MAD, phase/indexing, grain ID, and a GB overlay. The overlay uses visually distinct 2-5, 5-15, and at least 15 degree boundaries. Add HCP extension-twin candidate overlays only as angular-axis candidates and label them as non-unique.

**Step 5: Implement boundary network tables**

Export every eligible segment with physical length, angle, axis when stable, source class, and threshold bin. Summaries report total line density, LAGB/HAGB fractions by length, number-weighted sensitivity, cumulative distributions, and detection-floor sensitivity at 0.5, 1, 2, and 5 degrees.

**Step 6: Run full-data test**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath('C:/Users/22069/Documents/MATLAB/mtex-6.1.1'); startup_mtex; addpath('tools/mtex'); test_comprehensive_maps_morphology_boundaries('data/ebsd_kpl_250221_7_df/scans','.codex_tmp/ebsd_maps_test')"
```

Expected: no nonlocal native-grid boundaries, boundary-length conservation within numerical tolerance, no raw-data writes, and all artifacts present.

### Task 4: Quantify intragranular deformation and texture

**Files:**

- Create: `tools/mtex/compute_intragranular_metrics.m`
- Create: `tools/mtex/compute_texture_metrics.m`
- Create: `tools/mtex/generate_comprehensive_intragranular_texture.m`
- Create: `tools/mtex/test_comprehensive_intragranular_texture.m`
- Reuse: `tools/mtex/compute_c_axis_ad_statistics.m`
- Reuse: `tools/mtex/c_axis_angles_to_ad.m`

**Step 1: Write invariant tests**

On a constant-orientation synthetic grid, assert KAM = 0, GROD angle = 0, GOS = 0, c-axis-to-AD angle is bounded by 0-90 degrees, and an approximately uniform orientation set has lower texture concentration than a single-orientation set.

**Step 2: Run and observe failure**

Use the same MTEX startup command as Task 3 and call `test_comprehensive_intragranular_texture`.

**Step 3: Implement intragranular quantities**

Use native-grid `KAM` with grain IDs, primary order 1 and thresholds 2 and 5 degrees; add order 2 sensitivity. Compute `calcGROD(ebsd,grains)` angle and specimen/crystal-axis summaries where numerically stable, plus grain `GOS`. Export per-pixel tables only for registered metrics needed for reproducibility; avoid duplicating raw Euler data.

**Step 4: Implement texture quantities**

Calculate area-weighted grain-mean and pixel-weighted ODFs separately. Export {0001}, {10-10}, and {11-20} pole figures; AD/TD/ND inverse pole figures; consistent ODF sections; maximum MRD, `norm(odf)^2` texture index, `calcMIndex(odf)`, and `entropy(odf)`. Export c-axis acute angle to AD, full-sphere components, and azimuth about AD. Report raw-to-denoised ODF distances using a fixed common grid/kernel.

**Step 5: Run full-data test**

Expected: finite summary values, identical plotting conventions across deformation states, raw and denoised distinctly labeled, and all registered texture metrics reproducible from exported parameters.

### Task 5: Compute axial slip/twin propensity

**Files:**

- Create: `tools/mtex/compute_axial_propensity_metrics.m`
- Create: `tools/mtex/generate_comprehensive_axial_propensity.m`
- Create: `tools/mtex/test_comprehensive_axial_propensity.m`

**Step 1: Write crystallographic bounds tests**

For known HCP orientations, assert absolute Schmid factors lie in 0-0.5 and that symmetrising a family cannot reduce its maximum Schmid factor. Test basal `<a>`, prismatic `<a>`, pyramidal `<a>`, pyramidal `<c+a>`, and explicitly defined extension/contraction twin systems if their MTEX construction is validated.

**Step 2: Run and observe failure**

Use the MTEX startup command and call `test_comprehensive_axial_propensity`.

**Step 3: Implement grainwise Schmid summaries**

Rotate each HCP slip/twin family to specimen coordinates and evaluate uniaxial tension parallel to `vector3d.X`. Export the maximum absolute Schmid factor per family, winning variant index, area weight, and orientation. Do not assign CRSS-dependent activity unless a cited CRSS set is explicitly passed.

**Step 4: Add optional Taylor-factor sensitivity**

Use `calcTaylor(inv(ori)*eps,sS.symmetrise)` only for clearly documented slip-family sets and assumed equal or supplied CRSS values. Export this as sensitivity, not as a primary mechanism discriminator.

**Step 5: Run full-data test**

Expected: bound checks pass, every retained grain is represented exactly once per family, area weights sum correctly, and raw/denoised comparisons are paired by sample.

### Task 6: Integrate modules with tensile data and produce the manuscript-facing bundle

**Files:**

- Create: `tools/mtex/generate_comprehensive_raw_denoised_comparison.m`
- Create: `tools/mtex/generate_comprehensive_tensile_integration.m`
- Create: `tools/mtex/run_comprehensive_ebsd_analysis.m`
- Create: `tools/mtex/test_comprehensive_ebsd_analysis.m`
- Create: `results/mtex_ebsd_comprehensive/README.md` through the runner
- Create: `results/mtex_ebsd_comprehensive/analysis_manifest.csv` through the runner

**Step 1: Write integration tests**

Assert exact joins by diameter/cold reduction, preserve tensile `yield_status`, never average rows flagged `not_reliable_no_elastic_segment_or_offset_crossing` into Rp0.2 summaries, and report replicate count and dispersion for each tensile quantity.

**Step 2: Implement raw/denoised comparison**

Merge every paired EBSD summary metric, export absolute and relative differences, rank agreement, and direction-of-trend agreement. Flag metrics whose scientific interpretation changes under denoising.

**Step 3: Implement tensile integration**

Read `gr4b23271_cold_deformation_tensile_summary.csv` and the user-exclusion table. Aggregate UTS, reliable Rp0.2, uniform elongation, and true-stress/true-strain work-hardening descriptors with replicate uncertainty. With only six EBSD states, report Spearman correlations and leave-one-state-out sensitivity as descriptive associations; do not fit a high-parameter predictive model.

**Step 4: Implement the top-level runner**

```matlab
run_comprehensive_ebsd_analysis(projectRoot, outputRoot)
```

The runner creates a fresh derived-output tree, calls Tasks 2-5 in dependency order, verifies the output contract, writes exact input SHA-256 hashes and analysis parameters, and records the MATLAB/MTEX versions. It never deletes outside the explicit `outputRoot`.

**Step 5: Generate the README**

The README explains coordinate definitions, registered primary metrics, sensitivity analyses, raw/denoised roles, exclusions, output inventory, and interpretation limits in manuscript-ready Chinese.

**Step 6: Run integration test**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath('C:/Users/22069/Documents/MATLAB/mtex-6.1.1'); startup_mtex; addpath('tools/mtex'); test_comprehensive_ebsd_analysis(pwd,'.codex_tmp/ebsd_full_test')"
```

Expected: every contract artifact exists, table keys are unique, joins are complete, and no source hash changes.

### Task 7: Full production run, scientific review, and commit

**Step 1: Record input hashes before execution**

```powershell
Get-ChildItem 'data/ebsd_kpl_250221_7_df/scans' -Recurse -Filter '*.ctf' | Get-FileHash -Algorithm SHA256 | Sort-Object Path | Export-Csv '.codex_tmp/ctf_hashes_before.csv' -NoTypeInformation
```

**Step 2: Run all fast tests**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath('C:/Users/22069/Documents/MATLAB/mtex-6.1.1'); startup_mtex; addpath('tools/mtex'); test_comprehensive_ebsd_contract; test_comprehensive_ebsd_foundation; test_comprehensive_maps_morphology_boundaries; test_comprehensive_intragranular_texture; test_comprehensive_axial_propensity"
```

Expected: all tests print `passed` and MATLAB exits 0.

**Step 3: Run the production analysis**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath('C:/Users/22069/Documents/MATLAB/mtex-6.1.1'); startup_mtex; addpath('tools/mtex'); run_comprehensive_ebsd_analysis(pwd,fullfile(pwd,'results','mtex_ebsd_comprehensive'))"
```

**Step 4: Verify raw-data preservation and artifact completeness**

Recompute hashes and compare to the before file. Run the full integration test against the production output. Inspect every composite PNG for consistent orientation, readable legends, non-clipped labels, and common color limits.

**Step 5: Perform the scientific interpretation review**

Create a short derived note in `results/mtex_ebsd_comprehensive/README.md` that classifies each mechanism statement as directly measured, orientation-derived, proxy-based, or literature-dependent. Explicitly test whether the data support: early c-axis redistribution, later subboundary/intragranular misorientation increase, both concurrently, or neither monotonically.

**Step 6: Commit independently implemented packages and final integration**

```powershell
git add tools/mtex results/mtex_ebsd_comprehensive
git commit -m "feat: add comprehensive EBSD analysis workflow"
```

## Parallel agent ownership after Task 2

- Agent A owns Task 3 files only.
- Agent B owns Task 4 files only.
- Agent C owns Task 5 files only.
- The root agent owns Tasks 1, 2, 6, and 7, reviews all agent changes, resolves shared-interface issues, and performs production verification.
- Parallel agents must not modify shared foundation files, the output contract, tensile inputs, or raw CTF data. They must not commit while other agents are editing the shared worktree; the root agent stages and commits reviewed changes.

## Completion criteria

- All 12 CTF inputs are auditable and all six raw/denoised pairs are spatially identical.
- Every quantitative claim can be traced to an exported CSV and registered parameter set.
- Raw data remain byte-identical.
- Original-data trends and denoising sensitivity are separately visible.
- The horizontal map axis is consistently labeled AD in code, tables, figures, and prose.
- The final interpretation distinguishes direct observations from mechanism inference and does not force the supervisor's proposed two-stage mechanism when the measured trends do not support it.
