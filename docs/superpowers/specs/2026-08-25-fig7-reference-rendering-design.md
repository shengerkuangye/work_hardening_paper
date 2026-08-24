# Fig. 7 Reference-Style Texture Rendering Design

## Objective

Create a second rendering of the six Gr4B23271 EBSD texture states whose visual effect closely follows Fig. 7 of Xia et al. The change concerns plotting and composition only. The input CTF files, Ti-Hex phase selection, ODF calculation, kernel halfwidth, sampling resolution, and coordinate convention remain unchanged.

The existing common-scale contour figures remain available as the quantitative comparison version. The new outputs use distinct filenames and do not overwrite that version.

## Reference characteristics to reproduce

- Three circular pole figures in one row: `{0001}`, `{10-10}`, and `{11-20}`.
- Smooth continuous density fields without black contour lines.
- A blue-green-yellow-red density palette close to the reference figure.
- Black circular outlines and specimen-coordinate crosshairs in every pole figure.
- `AD` at the upper end of the vertical crosshair and `TD/RD` at the left end of the horizontal crosshair. These labels preserve the registered coordinate convention of the present EBSD data rather than copying an inapplicable direction name.
- Three inverse pole figures in the second row for `AD`, `TD/RD`, and `ND`.
- One vertical color scale at the right of each row with explicit `Max=...` and `Min=0.00` labels.
- Compact serif-style labels, reduced whitespace, and a thin dashed grouping frame similar to the published composition.

## Scale policy

The reference-style figures use the actual maximum density of each state's PF row and IPF row, matching the visual behavior of the reference article. This makes peaks and spatial distributions easy to see but prevents direct comparison by color alone between different states.

The existing common 0–7 MRD figures are therefore retained for cross-state quantitative comparison. The summary CSV records both state-specific maxima and the common-scale maxima.

## Output structure

For each of the six deformation states, generate one compact reference-style PNG containing the PF and IPF rows. Generate an additional 2 × 3 montage in PNG, TIFF, and PDF formats. The montage uses panel labels `(a)`–`(f)` for the six deformation states while keeping each state visually self-contained.

All outputs are written under a new `reference_style` subdirectory of `results/mtex_fig7_texture_gallery` and exported at 600 dpi.

## Verification

- Confirm six state figures and one montage in all required formats.
- Confirm all PNG files report 600 dpi.
- Confirm no black contour lines are present inside the density fields.
- Confirm PF crosshairs and direction labels are visible but do not cover density maxima.
- Confirm every PF/IPF row displays its own numeric maximum and `Min=0.00`.
- Visually compare the undeformed and maximum-deformation state figures with the published Fig. 7.
- Confirm MATLAB completes without stderr output or residual run-owned processes.

