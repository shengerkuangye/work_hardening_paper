# MTEX 6.1.1 Installation and CTF Smoke-Test Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Install the official MTEX 6.1.1 Windows release for MATLAB R2025a, provide a project-local launcher, and verify that the TA4-M CTF can produce a temporary α-Ti `{0001}` pole figure without modifying raw data.

**Architecture:** Keep the third-party MTEX distribution outside the Git repository at `C:\Users\22069\Documents\MATLAB\mtex-6.1.1`. Keep only a focused project launcher and a repeatable CTF smoke-test function under `tools/mtex/`; invoke MATLAB by absolute path so the workflow does not depend on the system PATH.

**Tech Stack:** Windows PowerShell 5.1, MATLAB R2025a, MTEX 6.1.1, Oxford/Channel 5 CTF.

## Global Constraints

- Follow `docs/superpowers/specs/2026-07-19-mtex-alpha-ti-c-axis-pole-figure-design.md`.
- Install exactly MTEX 6.1.1 at `C:\Users\22069\Documents\MATLAB\mtex-6.1.1`.
- Download the official binary release asset, not GitHub's source archive.
- Verify SHA-256 `5a9d0bcd771db0ef88125f3e368f1848e1aca3a9734efeb50593624aa31cedad` before extraction.
- Invoke MATLAB using `C:\Program Files\MATLAB\R2025a\bin\matlab.exe`.
- Do not overwrite or export over any `.ctf` file.
- Use the original TA4-M CTF for the smoke test; exclude unindexed and cubic-Ti pixels from the temporary ODF.
- Keep MTEX itself outside Git; commit only project scripts and documentation.

---

## File Structure

- External install: `C:\Users\22069\Documents\MATLAB\mtex-6.1.1\` — official MTEX distribution.
- Create: `tools/mtex/startup_project_mtex.m` — deterministic project-local MTEX launcher.
- Create: `tools/mtex/verify_mtex_ctf.m` — imports TA4-M CTF, validates phase/grid metadata, and writes a temporary `{0001}` pole figure.
- Temporary output: `.codex_tmp/mtex-smoke/ta4_m_0001_pf.png` — smoke-test artifact, not a manuscript result.

### Task 1: Install and verify the official MTEX release

**Files:**
- Create externally: `C:\Users\22069\Documents\MATLAB\mtex-6.1.1\`
- Create temporarily: `C:\Users\22069\AppData\Local\Temp\codex-mtex-6.1.1\mtex-6.1.1.zip`

**Interfaces:**
- Consumes: Official release asset URL and SHA-256 from the Global Constraints.
- Produces: A complete MTEX tree whose `VERSION` file reports `6.1.1` and whose `startup_mtex.m` is present.

- [ ] **Step 1: Run the installation preflight**

```powershell
$install = 'C:\Users\22069\Documents\MATLAB\mtex-6.1.1'
$matlab = 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe'
[pscustomobject]@{
  InstallAlreadyExists = Test-Path -LiteralPath $install
  MatlabExists = Test-Path -LiteralPath $matlab
  FreeGB = [math]::Round((Get-PSDrive C).Free / 1GB, 1)
}
```

Expected: `InstallAlreadyExists=False`, `MatlabExists=True`, and `FreeGB` comfortably exceeds 1 GB. If the install directory already exists, stop and inspect it; do not merge or overwrite an unknown installation.

- [ ] **Step 2: Download the official binary ZIP**

```powershell
$tempRoot = 'C:\Users\22069\AppData\Local\Temp\codex-mtex-6.1.1'
$zip = Join-Path $tempRoot 'mtex-6.1.1.zip'
New-Item -ItemType Directory -Force -Path $tempRoot | Out-Null
Invoke-WebRequest -Uri 'https://github.com/mtex-toolbox/mtex/releases/download/mtex-6.1.1/mtex-6.1.1.zip' -OutFile $zip
Get-Item -LiteralPath $zip | Select-Object FullName,Length
```

Expected: file length `196424572` bytes.

- [ ] **Step 3: Verify the release digest before extraction**

```powershell
$actual = (Get-FileHash -Algorithm SHA256 -LiteralPath $zip).Hash.ToLowerInvariant()
$expected = '5a9d0bcd771db0ef88125f3e368f1848e1aca3a9734efeb50593624aa31cedad'
if ($actual -ne $expected) { throw "MTEX checksum mismatch: $actual" }
$actual
```

Expected: exactly `5a9d0bcd771db0ef88125f3e368f1848e1aca3a9734efeb50593624aa31cedad`.

- [ ] **Step 4: Inspect the archive root and extract without overwriting**

```powershell
tar -tf $zip | Select-Object -First 10
$matlabTools = 'C:\Users\22069\Documents\MATLAB'
New-Item -ItemType Directory -Force -Path $matlabTools | Out-Null
if (Test-Path -LiteralPath $install) { throw "Refusing to overwrite $install" }
Expand-Archive -LiteralPath $zip -DestinationPath $matlabTools
```

Expected: the archive contains the `mtex-6.1.1/` root and extraction creates the exact installation directory.

- [ ] **Step 5: Verify essential installation files**

```powershell
$required = @(
  (Join-Path $install 'startup_mtex.m'),
  (Join-Path $install 'VERSION'),
  (Join-Path $install 'EBSDAnalysis'),
  (Join-Path $install 'interfaces'),
  (Join-Path $install 'mex')
)
$required | ForEach-Object { if (-not (Test-Path -LiteralPath $_)) { throw "Missing MTEX component: $_" } }
(Get-Content -Raw -LiteralPath (Join-Path $install 'VERSION')).Trim()
```

Expected: version output `6.1.1`.

### Task 2: Add and test a project-local MTEX launcher

**Files:**
- Create: `tools/mtex/startup_project_mtex.m`

**Interfaces:**
- Consumes: External MTEX root `C:\Users\22069\Documents\MATLAB\mtex-6.1.1`.
- Produces: `mtexRoot = startup_project_mtex()` returning a MATLAB string scalar and asserting that MTEX resolves from that root.

- [ ] **Step 1: Run the launcher test before creating the function**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath(fullfile(pwd,'tools','mtex')); startup_project_mtex"
```

Expected: MATLAB fails with `Unrecognized function or variable 'startup_project_mtex'`.

- [ ] **Step 2: Create the minimal deterministic launcher**

```matlab
function mtexRoot = startup_project_mtex()
%STARTUP_PROJECT_MTEX Load the project-pinned MTEX installation.

mtexRoot = "C:\Users\22069\Documents\MATLAB\mtex-6.1.1";
assert(isfolder(mtexRoot), "MTEX installation not found: %s", mtexRoot);

oldFolder = pwd;
restoreFolder = onCleanup(@() cd(oldFolder)); %#ok<NASGU>
cd(mtexRoot);
startup_mtex;

ebsdPath = string(which("EBSD"));
densityPath = string(which("calcDensity"));
assert(strlength(ebsdPath) > 0, "MTEX EBSD class was not loaded.");
assert(strlength(densityPath) > 0, "MTEX calcDensity was not loaded.");
assert(startsWith(lower(ebsdPath), lower(mtexRoot)), ...
  "EBSD resolves outside the pinned MTEX root: %s", ebsdPath);
assert(startsWith(lower(densityPath), lower(mtexRoot)), ...
  "calcDensity resolves outside the pinned MTEX root: %s", densityPath);
end
```

- [ ] **Step 3: Run the launcher test**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath(fullfile(pwd,'tools','mtex')); root=startup_project_mtex; fprintf('MTEX_ROOT=%s\n',root); fprintf('MTEX_VERSION=%s\n',mtex_version); fprintf('EBSD_PATH=%s\n',which('EBSD')); fprintf('CALCDENSITY_PATH=%s\n',which('calcDensity'));"
```

Expected: exit code 0, `MTEX_VERSION=6.1.1`, and both resolved paths begin with the pinned MTEX root.

- [ ] **Step 4: Commit the launcher**

```powershell
git add -- tools/mtex/startup_project_mtex.m
git commit -m "build: add project MTEX launcher"
```

### Task 3: Add a CTF and `{0001}` pole-figure smoke test

**Files:**
- Create: `tools/mtex/verify_mtex_ctf.m`
- Create temporarily at runtime: `.codex_tmp/mtex-smoke/ta4_m_0001_pf.png`

**Interfaces:**
- Consumes: `startup_project_mtex()` and the original TA4-M CTF.
- Produces: `report = verify_mtex_ctf(projectRoot)` with fields `mtexVersion`, `pixelCount`, `alphaPixelCount`, `scanStepUm`, `phaseNames`, and `outputFile`.

- [ ] **Step 1: Record the raw CTF hash and run the smoke test before implementation**

```powershell
$ctf = 'data\ebsd_kpl_250221_7_df\scans\d7\ebsd_sample_7_map_15.ctf'
Get-FileHash -Algorithm SHA256 -LiteralPath $ctf
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath(fullfile(pwd,'tools','mtex')); verify_mtex_ctf(pwd)"
```

Expected: MATLAB fails with `Unrecognized function or variable 'verify_mtex_ctf'`; retain the printed SHA-256 for Step 4.

- [ ] **Step 2: Create the CTF smoke-test function**

```matlab
function report = verify_mtex_ctf(projectRoot)
%VERIFY_MTEX_CTF Validate MTEX by importing TA4-M and plotting {0001}.

arguments
  projectRoot (1,1) string = string(pwd)
end

startup_project_mtex();

ctfFile = fullfile(projectRoot, "data", "ebsd_kpl_250221_7_df", ...
  "scans", "d7", "ebsd_sample_7_map_15.ctf");
assert(isfile(ctfFile), "CTF file not found: %s", ctfFile);

ebsd = EBSD.load(ctfFile, "convertEuler2SpatialReferenceFrame");
assert(length(ebsd) == 600 * 600, ...
  "Unexpected pixel count: %d", length(ebsd));

alpha = ebsd("Ti-Hex");
assert(~isempty(alpha), "Ti-Hex phase was not imported.");
assert(length(alpha) > 350000, ...
  "Unexpected Ti-Hex pixel count: %d", length(alpha));

xValues = unique(ebsd.x);
yValues = unique(ebsd.y);
dx = median(diff(xValues));
dy = median(diff(yValues));
assert(abs(dx - 0.5) < 1e-10 && abs(dy - 0.5) < 1e-10, ...
  "Unexpected scan step: dx=%g, dy=%g", dx, dy);

psi = deLaValleePoussinKernel("halfwidth", 5 * degree);
odf = calcDensity(alpha.orientations, "kernel", psi);
cAxis = Miller(0, 0, 0, 1, alpha.CS);

outputDir = fullfile(projectRoot, ".codex_tmp", "mtex-smoke");
if ~isfolder(outputDir), mkdir(outputDir); end
outputFile = fullfile(outputDir, "ta4_m_0001_pf.png");

fig = figure("Visible", "off");
cleanupFigure = onCleanup(@() close(fig)); %#ok<NASGU>
plotPDF(odf, cAxis, "antipodal", "contourf");
mtexColorbar("title", "m.r.d.");
saveFigure(outputFile);
assert(isfile(outputFile), "Pole-figure output was not created.");

phaseNames = strings(numel(ebsd.CSList), 1);
for k = 1:numel(ebsd.CSList)
  phaseNames(k) = string(ebsd.CSList{k}.mineral);
end

report = struct( ...
  "mtexVersion", string(mtex_version), ...
  "pixelCount", length(ebsd), ...
  "alphaPixelCount", length(alpha), ...
  "scanStepUm", [dx, dy], ...
  "phaseNames", phaseNames, ...
  "outputFile", string(outputFile));

fprintf("MTEX_VERSION=%s\n", report.mtexVersion);
fprintf("PIXELS=%d\n", report.pixelCount);
fprintf("ALPHA_PIXELS=%d\n", report.alphaPixelCount);
fprintf("SCAN_STEP_UM=%.4f,%.4f\n", report.scanStepUm);
fprintf("OUTPUT=%s\n", report.outputFile);
end
```

- [ ] **Step 3: Run the CTF smoke test**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath(fullfile(pwd,'tools','mtex')); r=verify_mtex_ctf(string(pwd)); assert(r.pixelCount==360000); assert(r.alphaPixelCount==358203); assert(all(abs(r.scanStepUm-[0.5 0.5])<1e-10));"
```

Expected: exit code 0, `PIXELS=360000`, `ALPHA_PIXELS=358203`, `SCAN_STEP_UM=0.5000,0.5000`, and a PNG path under `.codex_tmp/mtex-smoke/`.

- [ ] **Step 4: Verify the raw CTF is unchanged and inspect the output**

```powershell
Get-FileHash -Algorithm SHA256 -LiteralPath $ctf
Get-Item -LiteralPath '.codex_tmp\mtex-smoke\ta4_m_0001_pf.png' | Select-Object FullName,Length
```

Expected: the CTF SHA-256 exactly matches Step 1 and the PNG has nonzero length. Open the PNG for visual inspection; the plot must be legible and free of MATLAB error annotations.

- [ ] **Step 5: Commit the smoke-test function**

```powershell
git add -- tools/mtex/verify_mtex_ctf.m
git commit -m "test: verify MTEX CTF import and pole figure"
```

### Task 4: Final installation verification

**Files:**
- Verify: `tools/mtex/startup_project_mtex.m`
- Verify: `tools/mtex/verify_mtex_ctf.m`
- Verify externally: `C:\Users\22069\Documents\MATLAB\mtex-6.1.1\`

**Interfaces:**
- Consumes: Tasks 1–3.
- Produces: Evidence that MATLAB, MTEX, the project launcher, CTF import, α-Ti selection, and temporary `{0001}` plotting all work together.

- [ ] **Step 1: Run the complete verification from a fresh MATLAB process**

```powershell
& 'C:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "restoredefaultpath; addpath(fullfile(pwd,'tools','mtex')); root=startup_project_mtex; r=verify_mtex_ctf(string(pwd)); assert(string(mtex_version)=='6.1.1'); assert(startsWith(lower(string(which('EBSD'))),lower(root))); assert(r.pixelCount==360000);"
```

Expected: exit code 0 with no missing MEX, path-conflict, CTF-import, or plotting errors.

- [ ] **Step 2: Confirm repository state and installation boundary**

```powershell
git status --short
git log -3 --oneline
Get-ChildItem -LiteralPath 'C:\Users\22069\Documents\MATLAB\mtex-6.1.1' | Select-Object -First 10 Name
```

Expected: only intentional project files are tracked, MTEX itself is outside the repository, and the launcher and smoke-test commits are present.
