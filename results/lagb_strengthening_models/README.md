# LAGB strengthening models

## Model scope

The calculation follows the literature decomposition `sigma_y = sigma_0 + sigma_HAGB + sigma_LAGB`. KAM-derived GND strengthening is not added separately because the LAGB model already represents dislocations arranged in low-angle boundaries.

Strict model: literature-form equation and independent CP-alpha-Ti parameters; the 7 mm condition calibrates the matrix/base term only.

Closest model: the same structure, with one dimensionless multiplier `q` on the LAGB term chosen to minimize MAPE over the five deformed conditions.

## Parameters and metrics

- Strict baseline: 405.721 MPa
- Closest q: 1.871509
- Closest baseline: 374.456 MPa
- Strict MAPE, deformed only: 18.673%
- Closest MAPE, deformed only: 3.886%
- Closest leave-one-deformed-state-out MAPE: 4.161%

## Results

| Diameter (mm) | Experimental (MPa) | HAGB (MPa) | LAGB strict (MPa) | Strict prediction (MPa) | Strict deviation (%) | Closest prediction (MPa) | Closest deviation (%) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 7.00 | 580.5 | 138.9 | 35.9 | 580.5 | +0.00 | 580.5 | +0.00 |
| 6.48 | 845.0 | 132.0 | 197.0 | 734.7 | -13.05 | 875.1 | +3.56 |
| 6.02 | 900.0 | 133.2 | 140.8 | 679.8 | -24.47 | 771.2 | -14.31 |
| 5.60 | 927.0 | 132.9 | 224.2 | 762.9 | -17.71 | 927.0 | +0.00 |
| 5.25 | 982.5 | 138.7 | 254.6 | 799.1 | -18.67 | 989.7 | +0.73 |
| 5.00 | 980.0 | 138.0 | 245.5 | 789.2 | -19.47 | 971.9 | -0.82 |

## Interpretation boundary

The strict version is the defensible literature-form estimate with currently available EBSD quantities. Its stereological intercept assumes a two-dimensional isotropic boundary network: `L = pi/(2 lambda_LAGB)`. Direct directional line-intercept measurements should replace this approximation for the elongated cold-worked structure.

The closest version is a fitted effective model, not independent validation. The multiplier can represent boundary anisotropy, line-intercept conversion, alpha, Taylor-factor, and unresolved substructure effects; it must not be assigned uniquely to any one mechanism without additional measurements.

Sources:

- https://doi.org/10.2320/matertrans.M2018033
- https://doi.org/10.1016/S1359-6454(97)00365-0
