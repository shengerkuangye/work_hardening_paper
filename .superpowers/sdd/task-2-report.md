# Task 2 Report: Seven-section alpha-Ti ODF diagnostics

## Scope and result

Implemented `tools/mtex/generate_odf_diameter_montage.m` for the six raw
Ti-Hex EBSD scans. The function now:

- selects the six raw catalog rows in the required diameter order;
- uses each imported orientation array's existing CS/SS;
- validates `6/mmm`, rotational proper group `622`, specimen symmetry `1`,
  and the MTEX `fundamentalRegionEuler` limits `360/90/60 deg`;
- calculates an equal-pixel-weighted ODF with a 5 deg
  `SO3DeLaValleePoussinKernel` and normalizes its mean to one;
- converts the ODF to `SO3FunHarmonic` and renormalizes its mean to one;
- evaluates every ODF on the literal
  `phi1 = 0:5:355`, `Phi = 0:5:90`, and `phi2 = 0:10:60 deg` grids;
- returns 6-row and 42-row tables and writes the two required CSV files.

PNG/PDF rendering was intentionally not implemented in this task.

## TDD and verification evidence

All commands were run from:

`C:\Users\Admin\Desktop\GithubRepo\work_hardening_paper-odf-sections`

### RED baseline

Command:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "startup_mtex('noMenu'); addpath('tools/mtex'); test_generate_odf_diameter_montage(string(fullfile(pwd,'data','ebsd_kpl_250221_7_df','scans')), string(fullfile(pwd,'.codex_tmp','odf-montage-summary-green')));"
```

- Wall time: 12.83 s
- Exit: 1
- Expected first failure:

```text
函数或变量 'generate_odf_diameter_montage' 无法识别。
出错 test_generate_odf_diameter_montage (第 5 行)
```

This established that the new entry point was absent before implementation.

### First calculation run

The same test command was run after the initial implementation.

- Wall time: 37.48 s
- Exit: 1
- All assertions through the two CSV-existence checks passed.
- The first failure was the intentionally unimplemented PNG assertion:

```text
出错 test_generate_odf_diameter_montage (第 55 行)
assert(isfile(fullfile(outputRoot,"odf_diameter_full_sections.png")));
```

The cold-reduction expression was then simplified without changing behavior,
so a fresh final run was required.

### Final fresh calculation run

The same full test command was run against the final source.

- Wall time: 48.99 s
- Exit: 1, expected for Task 2
- The first failure remained exactly the line-55 PNG assertion shown above.
- Therefore the preceding symmetry, ordering, MRD, global color-limit, table,
  and CSV-existence assertions all executed and passed.
- The PDF assertion was not reached because MATLAB stopped at the earlier PNG
  assertion. Neither image file exists, as required for this task boundary.

### CSV round-trip validation

The test's round-trip assertions occur after its image assertions and are
therefore unreachable in the expected Task 2 failure. They were executed
separately against the CSV files from the final fresh calculation run:

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "outputRoot=string(fullfile(pwd,'.codex_tmp','odf-montage-summary-green')); s=readtable(fullfile(outputRoot,'odf_diameter_summary.csv'),'TextType','string'); q=readtable(fullfile(outputRoot,'odf_diameter_section_summary.csv'),'TextType','string'); sampleColumns={'sample','diameter_mm','cold_reduction_percent','input_path','valid_ti_hex_orientation_count','crystal_symmetry','rotational_fundamental_zone','specimen_symmetry','kernel_halfwidth_deg','grid_resolution_deg','maximum_section_mrd','global_color_limit_max_mrd'}; sectionColumns={'sample','diameter_mm','cold_reduction_percent','input_path','valid_ti_hex_orientation_count','crystal_symmetry','rotational_fundamental_zone','specimen_symmetry','kernel_halfwidth_deg','grid_resolution_deg','phi2_deg','section_maximum_mrd','global_color_limit_max_mrd'}; assert(isequal(s.Properties.VariableNames,sampleColumns)); assert(isequal(q.Properties.VariableNames,sectionColumns)); assert(height(s)==6 && height(q)==42); assert(q.global_color_limit_max_mrd(1)==max(q.section_maximum_mrd)); assert(all(s.global_color_limit_max_mrd==q.global_color_limit_max_mrd(1))); disp('CSV_ROUNDTRIP_PASS');"
```

- Wall time: 15.64 s
- Exit: 0
- Output: `CSV_ROUNDTRIP_PASS`

An earlier auxiliary round-trip invocation exited 1 before any assertion
because shell quoting produced an incomplete MATLAB statement. It was
discarded and replaced by the successful command above.

## Numerical summary

The global seven-section color-limit maximum is
`10.4912388697768 MRD`.

| Sample | Diameter (mm) | Cold reduction (%) | Valid Ti-Hex orientations | Maximum of 7 sections (MRD) |
|---|---:|---:|---:|---:|
| 7d | 7.00 | 0 | 358203 | 9.89035261390408 |
| 6.48d | 6.48 | 14.305306122449 | 350441 | 8.71699374087264 |
| 6.02d | 6.02 | 26.04 | 355266 | 7.16900270789942 |
| 5.6d | 5.60 | 36 | 335389 | 7.74574847172371 |
| 5.25d | 5.25 | 43.75 | 317224 | 10.4912388697768 |
| 5d | 5.00 | 48.9795918367347 | 328099 | 8.52463479098233 |

The seven section maxima in `phi2 = 0,10,...,60 deg` order are:

- 7d: 7.719430, 9.890353, 7.116619, 6.143086, 6.714957,
  6.094915, 7.719430
- 6.48d: 8.716994, 7.929108, 7.465097, 7.235829, 7.026871,
  7.221306, 8.716994
- 6.02d: 7.169003, 6.648538, 6.847647, 7.114777, 7.160783,
  7.069651, 7.169003
- 5.6d: 7.583699, 7.745748, 7.521536, 6.682151, 7.157636,
  7.345393, 7.583699
- 5.25d: 10.491239, 10.105680, 9.128624, 9.244229, 9.092436,
  9.612214, 10.491239
- 5d: 8.524635, 8.447668, 6.798191, 6.798191, 6.798191,
  7.231588, 8.524635

## Self-review

- No local `fundamentalRegionEuler` function was defined. The implementation
  directly calls the MTEX function with three outputs.
- No new `crystalSymmetry` or `specimenSymmetry` object is instantiated, and
  no imported orientation symmetry is reassigned.
- `calcDensity` receives all imported Ti-Hex orientations with an explicit
  vector of unit weights.
- Both the RBF ODF and retained harmonic ODF are normalized to unit mean.
- The exact literal Euler grids and `real(eval(...))` evaluation path are
  used for all 42 sections.
- The sample and section table schemas, row counts, diameter-major ordering,
  and `phi2` ordering match the test contract.
- Cold reduction is calculated from area reduction relative to the 7 mm
  starting diameter because the shared catalog stores two-decimal rounded
  values for some samples.
- `git diff --check` reported no whitespace errors.
- The existing test file and all raw EBSD inputs were left unchanged.

## Remaining concern / intentional gap

The only intentional gap is PNG/PDF rendering, assigned to the subsequent
task. Consequently, the full test exits nonzero at the first image assertion;
all non-image checks, including separate CSV round-trip validation, pass.
