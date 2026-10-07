# Xia et al. (2026): TA16 cold rolling reference note

## Bibliographic record

Xia Y, Yu G, Yoon B-H, Cao P, Yi Y, Yu Z. Numerical simulation of cold working and microstructure-property regulation of TA16 titanium alloy. *Materials Science and Engineering: A*. 2026;975:150869. DOI: 10.1016/j.msea.2026.150869.

Source PDF: `references/xia_2026_ta16_cold_rolling_microstructure_property.pdf`.

## Material and process

- Material: TA16, nominal Ti-2Al-2.5Zr, not commercially pure TA4/Gr4B23271.
- Product/process: tube blank subjected to multi-pass reciprocating three-roll cold rolling.
- Methods: ABAQUS finite-element analysis, axial tensile testing, EBSD and TEM.
- EBSD step size: 0.7 μm; the paper reports OIM denoising before analysis.
- Simulation comparison: 25%, 35%, 45% and 55% reduction per pass; approximately 45% was selected as the preferred balance of stress uniformity and surface forming quality.

## Reported observations

- Mean grain size decreased from 19.68 μm in the initial material to 2.52 μm in the final state.
- The reported low-angle grain-boundary fraction reached 67.8%.
- A KAM-based expression was used to estimate GND density; the final reported value was 2.94 × 10^14 m^-2. This is an EBSD-derived estimate, not a direct measurement of total dislocation density.
- Texture evolved toward a strong {10-10}//AD fiber texture.
- The authors interpret the early stage as dominated by {10-12} tensile twinning and the later stage as dominated by dislocation slip, with basal and pyramidal <c+a> participation and suppressed prismatic slip under the evolved texture. TEM was used to supplement the EBSD/Schmid-factor interpretation.
- Final tensile strength was reported as 852 MPa and elongation as 9.5%.
- In the strengthening decomposition, dislocation strengthening was the largest contribution to total yield strength, whereas grain-refinement strengthening provided the larger increment between the compared intermediate and final rolled states.

## Use in the Gr4B23271 project

This paper is suitable for:

- organizing the process/load → tensile response → EBSD/TEM → deformation mechanism → strengthening contribution sequence;
- supporting the need to distinguish EBSD-derived orientation-gradient proxies from direct dislocation observations;
- illustrating why finite-element stress paths, texture, Schmid factors and TEM should be combined before assigning active slip or twin mechanisms;
- providing a comparison for the strength-ductility trade-off under severe cold working.

It must not be used as direct evidence that Gr4B23271/TA4 rotary-swaged bars undergo the same continuous grain refinement, twin transition, slip-system sequence or strengthening partition. The alloy chemistry, tube geometry, rolling kinematics, strain path, EBSD step size and available TEM evidence differ from the present project.
