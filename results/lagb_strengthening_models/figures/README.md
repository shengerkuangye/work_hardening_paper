# Fig. 15-style strengthening figures

Two figure families are generated from `lagb_strengthening_models.json`:

- `strict/individual/`: six state-specific figures for the literature-form model.
- `closest/individual/`: six state-specific figures for the one-parameter closest model.
- `strict_strengthening_montage.*`: six-panel strict-model montage.
- `closest_strengthening_montage.*`: six-panel closest-model montage.

Each calculated bar is stacked as matrix/base, HAGB strengthening and LAGB strengthening. The adjacent blue bar is the measured mean 0.2% proof stress. The predicted total and signed relative deviation are printed above the calculated bar.

SVG files are editable vector masters. The two montage PDF files are vector exports. PNG and LZW-compressed TIFF files are exported at 600 dpi. All panels use the same 0-1100 MPa y-axis range and the same colors.

The closest figure is a fitted effective model and must not be described as independent validation. The 6.02 mm state retains a -14.31% deviation.
