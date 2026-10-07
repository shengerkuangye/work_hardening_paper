# Reference library

## Newly registered project reference

- `xia_2026_ta16_cold_rolling_microstructure_property.pdf`
  - Xia Y, Yu G, Yoon B-H, Cao P, Yi Y, Yu Z.
  - *Numerical simulation of cold working and microstructure-property regulation of TA16 titanium alloy*.
  - Materials Science and Engineering: A 975 (2026) 150869.
  - DOI: 10.1016/j.msea.2026.150869.
  - Project role: comparison of the process-load-microstructure-property methodology and evidence hierarchy; TA16-specific mechanisms are not transferred directly to the Grade 4 CP-Ti specimens in this project.
  - Structured note: `extracted_literature_notes/xia_2026_ta16_cold_rolling_summary.md`.

## References moved out of temporary storage

- `acar_2017_ti7al_crystal_plasticity_odf.pdf` - ODF-based crystal-plasticity modeling of Ti-7Al; DOI 10.3390/met7110459.
- `wronski_2018_cp_ti_crystal_plasticity.pdf` - CP-Ti experiment and crystal-plasticity modeling; DOI 10.1016/j.msea.2018.03.017.
- `nemeth_2018_cp_ti_conform_ecap.pdf` - CP Ti Grade 2 microstructure, local texture and residual stress after CONFORM ECAP; DOI 10.3390/met8121000.
- `molodov_2007_cp_ti_annealing_texture.pdf` - cold-rolled CP-Ti annealing texture and grain growth; DOI 10.2320/matertrans.MI200701.
- `wang_2018_cp_ti_rotary_swaging_annealing.pdf` - CP Ti rotary swaging, annealing, texture and tensile/fatigue properties; DOI 10.3390/ma11112261.

`library.bib` registers the Xia et al. manuscript reference and the Xiong et al. CRSS calculation source. Other PDFs should be added to the BibTeX library only when they are cited in the text and their metadata have been checked.

## CRSS calculation source

- Xiong Y, Karamched P, Nguyen C-T, Collins DM, Magazzeni CM, Tarleton E, Wilkinson AJ. *Cold Creep of Titanium: Analysis of stress relaxation using synchrotron diffraction and crystal plasticity simulations*. Acta Materialia 199 (2020) 561–577. DOI: 10.1016/j.actamat.2020.08.010; arXiv:2003.01682.
- Project role: Table 4 provides the Grade 4 CP-Ti CRSS inputs used in `results/mtex_crss_normalized_slip_resistance/`: prismatic `<a>` 154 MPa, basal `<a>` 252 MPa and pyramidal `<c+a>` 756 MPa.
- Slip-family mapping: the paper discusses first- and second-order pyramidal slip and sets the pyramidal CRSS to three times the basal value; Table 4 does not independently calibrate the two orders. Applying 756 MPa to this project's MTEX first-order `pyramidalCA` geometry is an explicit comparison assumption. See [author manuscript, Sections 2.3/3.3 and Table 4](https://arxiv.org/pdf/2003.01682).
- Applicability limit: older project notes quoted O = 0.326–0.339 wt.% and Fe = 0.367–0.370 wt.%, but the composition certificate has not been traced to the tested bars. These ranges are not established specimen compositions and cannot substantiate a composition match to Xiong et al. The CRSS values remain literature model inputs for comparative screening, not material-specific measurements or calibration. See [material naming and traceability](/home/abcd/repos/work_hardening_paper/docs/material_naming_and_traceability.md).
