# Matrix + three-factor strengthening model

## Definition

The decomposition is `sigma_y = sigma_0 + sigma_SS + sigma_HP + sigma_dis`. The matrix friction stress is shown as an independent base layer; the three strengthening factors are solid solution, HAGB Hall--Petch, and KAM-GND dislocation strengthening.

## Parameters

- Matrix term: `172.5 MPa`.
- Solid-solution term: `258.11 MPa` from the O-equivalent Grade 4 CP-Ti model evaluation (`O_eq = 0.3484 wt.%`).
- HAGB term: `sigma_HP = 330/sqrt(d)` MPa, with `d` in micrometres.
- Dislocation term: `sigma_dis = alpha M G b sqrt(rho_GND)`, with `alpha = 0.2`, `M = 5.0`, `G = 45.6 GPa`, and `b = 0.295 nm`.

## Fit-free literature-parameter result

- Six-state MAPE: `13.12%`.
- Deformed-state MAPE: `13.47%`.
- Maximum absolute relative deviation: `16.81%`.
- RMSE: `124.8 MPa`.

## Bounded common-parameter calibration

To correct the two identified transfer mismatches without state-specific fitting, the calibrated version replaces the reference-alloy matrix term by one common effective matrix stress and converts KAM-GND to an effective total dislocation density through one common density ratio. Solid solution and HAGB parameters remain fixed.

- Calibrated matrix term: `76.75 MPa`.
- Dislocation strength multiplier: `1.854`.
- Effective `rho_total/rho_GND`: `3.437`.
- Six-state MAPE: `4.92%`.
- Deformed-state MAPE: `4.67%`.
- Maximum absolute relative deviation: `6.17%`.
- RMSE: `44.8 MPa`.

## Interpretation limits

- The solid-solution term is constant because all states come from the same material batch. It explains the alloy-composition baseline but not the increase with cold reduction.
- The O-equivalent literature model does not explicitly isolate the Fe contribution; therefore the numerical solid-solution term should be described as a literature-parameter estimate.
- The project O-equivalent value (0.3484 wt.%) is slightly above the source model's reported 0.14--0.32 wt.% validation interval. Its use is therefore a modest extrapolation, not an independently validated measurement of solid-solution strengthening.
- KAM supplies an EBSD-derived GND estimate rather than total dislocation density. The model does not add a separate LAGB term, avoiding direct double counting of boundary-arranged dislocations.
- The HAGB number-mean ECD changes little among the six states, so the Hall--Petch term is nearly constant.
- No parameter is fitted to the six experimental yield strengths. The remaining underprediction after deformation may include statistically stored dislocations, residual stress, texture-dependent Taylor response, and uncertainty in effective strengthening length scale; these are not assigned quantitatively here.
- The preceding statement applies to the literature-parameter version. The calibrated version uses two shared parameters fitted against all six yield strengths and is a descriptive closure, not an independent prediction or validation.

## Sources

- Solid-solution model: https://doi.org/10.1007/s12540-014-6004-8
- Independent CP-Ti Hall--Petch parameter: https://doi.org/10.2320/matertrans.M2018033
- Reference decomposition and matrix/Taylor parameters: attached Xia et al. reference paper.

## Figure outputs

- Six state figures: `figures/individual/` (editable SVG and 600 dpi PNG/TIFF).
- Six-state montage: `figures/reference_three_factor_strengthening_montage` (editable SVG, 600 dpi PNG/TIFF, and vector PDF).
- Calibrated six-state figures: `figures/calibrated/`, including individual figures and `reference_three_factor_strengthening_calibrated_montage` in SVG/PNG/TIFF/PDF.
