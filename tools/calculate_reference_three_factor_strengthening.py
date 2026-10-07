#!/usr/bin/env python3
"""Calculate the matrix + three-factor strengthening decomposition.

The three explicitly separated strengthening terms are solid solution,
HAGB Hall--Petch, and KAM-GND Taylor strengthening.  The matrix friction
stress is displayed separately and is not counted as a strengthening factor.
Raw experimental files are read only; all outputs are derived files.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TENSILE = ROOT / "data/tensile_data/gr4b23271_lab_tensile_by_diameter.csv"
HAGB = ROOT / "results/mtex_fig6_component_gallery/07_hagb_grain_size/hagb_grain_size_summary_used.csv"
GND = ROOT / "results/mtex_reference_style_fig14_gnd/reference_style_fig14_gnd_summary.csv"
OUTPUT = ROOT / "results/reference_three_factor_strengthening"

SAMPLES = ["7d", "6.48d", "6.02d", "5.6d", "5.25d", "5d"]
TENSILE_DIAMETERS = [7.0, 6.5, 6.02, 5.6, 5.25, 5.0]

# Literature-parameter model.  sigma_0 follows the reference paper's Fig. 15
# decomposition.  sigma_SS is the previously evaluated O_eq-based Grade 4
# CP-Ti value at 300 K and 1e-3 s^-1; it is kept separate from sigma_0.
SIGMA_MATRIX_MPA = 172.5
SIGMA_SOLID_SOLUTION_MPA = 258.11
O_EQ_WT_PERCENT = 0.3484
K_HAGB_MPA_UM_HALF = 330.0
TAYLOR_FACTOR = 5.0
ALPHA = 0.2
SHEAR_MODULUS_PA = 45.6e9
BURGERS_VECTOR_M = 0.295e-9


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def model_metrics(
    rows: list[dict[str, float | str]], prediction_key: str
) -> dict[str, float]:
    errors = [
        100.0
        * (float(row[prediction_key]) - float(row["yield_experimental_MPa"]))
        / float(row["yield_experimental_MPa"])
        for row in rows
    ]
    residuals = [
        float(row[prediction_key]) - float(row["yield_experimental_MPa"])
        for row in rows
    ]
    return {
        "mape_all_percent": sum(abs(value) for value in errors) / len(errors),
        "mape_deformed_only_percent": sum(abs(value) for value in errors[1:])
        / (len(errors) - 1),
        "max_absolute_relative_deviation_percent": max(abs(value) for value in errors),
        "rmse_MPa": math.sqrt(sum(value * value for value in residuals) / len(residuals)),
    }


def calibrate_common_parameters(
    rows: list[dict[str, float | str]],
) -> tuple[float, float]:
    """Minimize the largest relative deviation using two common parameters.

    The bounded grid keeps the true matrix term non-negative and the effective
    total-dislocation density between one and nine times the KAM-GND estimate.
    The fitted multiplier q acts on sigma_dis, so rho_total/rho_GND = q**2.
    No state receives an individual correction.
    """

    best: tuple[float, float, float] | None = None
    for matrix_step in range(0, 691):  # 0--172.5 MPa, 0.25 MPa step
        sigma_matrix = matrix_step * 0.25
        for multiplier_step in range(500, 1501):  # q=1.000--3.000
            q = multiplier_step * 0.002
            maximum_deviation = 0.0
            for row in rows:
                prediction = (
                    sigma_matrix
                    + float(row["sigma_solid_solution_MPa"])
                    + float(row["sigma_hagb_MPa"])
                    + q * float(row["sigma_dislocation_MPa"])
                )
                experimental = float(row["yield_experimental_MPa"])
                deviation = abs(100.0 * (prediction - experimental) / experimental)
                maximum_deviation = max(maximum_deviation, deviation)
            candidate = (maximum_deviation, sigma_matrix, q)
            if best is None or candidate < best:
                best = candidate
    assert best is not None
    return best[1], best[2]


def calculate() -> dict[str, object]:
    tensile = {
        float(row["nominal_diameter_group_mm"]): row for row in read_rows(TENSILE)
    }
    hagb = {row["sample"]: row for row in read_rows(HAGB)}
    gnd = {row["sample"]: row for row in read_rows(GND)}

    rows: list[dict[str, float | str]] = []
    for sample, tensile_diameter in zip(SAMPLES, TENSILE_DIAMETERS):
        h_row = hagb[sample]
        g_row = gnd[sample]
        t_row = tensile[tensile_diameter]
        grain_size_um = float(h_row["ecd_number_mean_um"])
        rho_gnd_1e14_m2 = float(g_row["gnd_mean_1e14_m2"])
        sigma_hagb = K_HAGB_MPA_UM_HALF / math.sqrt(grain_size_um)
        sigma_dislocation = (
            ALPHA
            * TAYLOR_FACTOR
            * SHEAR_MODULUS_PA
            * BURGERS_VECTOR_M
            * math.sqrt(rho_gnd_1e14_m2 * 1.0e14)
            / 1.0e6
        )
        prediction = (
            SIGMA_MATRIX_MPA
            + SIGMA_SOLID_SOLUTION_MPA
            + sigma_hagb
            + sigma_dislocation
        )
        experimental = float(t_row["Rp0.2_MPa_mean"])
        relative_deviation = 100.0 * (prediction - experimental) / experimental
        rows.append(
            {
                "sample": sample,
                "diameter_mm": float(h_row["diameter_mm"]),
                "cold_reduction_percent": float(h_row["cold_reduction_percent"]),
                "yield_experimental_MPa": experimental,
                "yield_std_MPa": float(t_row["Rp0.2_MPa_std"]),
                "hagb_number_mean_ecd_um": grain_size_um,
                "gnd_mean_1e14_m2": rho_gnd_1e14_m2,
                "sigma_matrix_MPa": SIGMA_MATRIX_MPA,
                "sigma_solid_solution_MPa": SIGMA_SOLID_SOLUTION_MPA,
                "sigma_hagb_MPa": sigma_hagb,
                "sigma_dislocation_MPa": sigma_dislocation,
                "yield_calculated_MPa": prediction,
                "relative_deviation_percent": relative_deviation,
            }
        )

    sigma_matrix_calibrated, dislocation_multiplier = calibrate_common_parameters(rows)
    density_ratio = dislocation_multiplier * dislocation_multiplier
    for row in rows:
        sigma_dislocation_calibrated = (
            dislocation_multiplier * float(row["sigma_dislocation_MPa"])
        )
        prediction_calibrated = (
            sigma_matrix_calibrated
            + float(row["sigma_solid_solution_MPa"])
            + float(row["sigma_hagb_MPa"])
            + sigma_dislocation_calibrated
        )
        experimental = float(row["yield_experimental_MPa"])
        row["sigma_matrix_calibrated_MPa"] = sigma_matrix_calibrated
        row["rho_total_to_gnd_ratio_calibrated"] = density_ratio
        row["sigma_dislocation_calibrated_MPa"] = sigma_dislocation_calibrated
        row["yield_calibrated_MPa"] = prediction_calibrated
        row["calibrated_relative_deviation_percent"] = (
            100.0 * (prediction_calibrated - experimental) / experimental
        )

    literature_metrics = model_metrics(rows, "yield_calculated_MPa")
    calibrated_metrics = model_metrics(rows, "yield_calibrated_MPa")

    return {
        "model_name": "matrix_plus_three_strengthening_factors",
        "equation": "sigma_y = sigma_0 + sigma_SS + sigma_HP + sigma_dis",
        "factor_definition": {
            "matrix": "sigma_0; displayed separately and not counted among the three strengthening factors",
            "factor_1": "O_eq-based solid-solution term sigma_SS",
            "factor_2": "HAGB Hall-Petch term sigma_HP",
            "factor_3": "KAM-GND Taylor term sigma_dis",
        },
        "parameters": {
            "sigma_matrix_MPa": SIGMA_MATRIX_MPA,
            "sigma_solid_solution_MPa": SIGMA_SOLID_SOLUTION_MPA,
            "oxygen_equivalent_wt_percent": O_EQ_WT_PERCENT,
            "k_hagb_MPa_um_half": K_HAGB_MPA_UM_HALF,
            "M": TAYLOR_FACTOR,
            "alpha": ALPHA,
            "G_GPa": SHEAR_MODULUS_PA / 1.0e9,
            "b_nm": BURGERS_VECTOR_M * 1.0e9,
        },
        "input_files": {
            "tensile": str(TENSILE.relative_to(ROOT)),
            "hagb_grain_size": str(HAGB.relative_to(ROOT)),
            "kam_gnd": str(GND.relative_to(ROOT)),
        },
        "source_urls": [
            "https://doi.org/10.1007/s12540-014-6004-8",
            "https://doi.org/10.2320/matertrans.M2018033",
        ],
        "calibration": {
            "name": "bounded_common_parameter_minimax",
            "objective": "minimum maximum absolute relative deviation across all six states",
            "state_specific_parameters": 0,
            "fitted_common_parameters": 2,
            "sigma_matrix_calibrated_MPa": sigma_matrix_calibrated,
            "dislocation_strength_multiplier": dislocation_multiplier,
            "rho_total_to_gnd_ratio": density_ratio,
            "equivalent_M_if_expressed_as_Taylor_factor": TAYLOR_FACTOR
            * dislocation_multiplier,
            "bounds": {
                "sigma_matrix_MPa": [0.0, SIGMA_MATRIX_MPA],
                "rho_total_to_gnd_ratio": [1.0, 9.0],
            },
            "grid_resolution": {
                "sigma_matrix_MPa": 0.25,
                "dislocation_strength_multiplier": 0.002,
            },
        },
        "literature_metrics": literature_metrics,
        "calibrated_metrics": calibrated_metrics,
        "rows": rows,
    }


def write_outputs(result: dict[str, object]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "reference_three_factor_strengthening.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    rows = result["rows"]
    csv_path = OUTPUT / "reference_three_factor_strengthening.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    literature_metrics = result["literature_metrics"]
    calibrated_metrics = result["calibrated_metrics"]
    calibration = result["calibration"]
    report = [
        "# Matrix + three-factor strengthening model",
        "",
        "## Definition",
        "",
        "The decomposition is `sigma_y = sigma_0 + sigma_SS + sigma_HP + sigma_dis`. "
        "The matrix friction stress is shown as an independent base layer; the three "
        "strengthening factors are solid solution, HAGB Hall--Petch, and KAM-GND "
        "dislocation strengthening.",
        "",
        "## Parameters",
        "",
        f"- Matrix term: `{SIGMA_MATRIX_MPA:.1f} MPa`.",
        f"- Solid-solution term: `{SIGMA_SOLID_SOLUTION_MPA:.2f} MPa` from the "
        f"O-equivalent Grade 4 CP-Ti model evaluation (`O_eq = {O_EQ_WT_PERCENT:.4f} wt.%`).",
        f"- HAGB term: `sigma_HP = {K_HAGB_MPA_UM_HALF:.0f}/sqrt(d)` MPa, with `d` in micrometres.",
        f"- Dislocation term: `sigma_dis = alpha M G b sqrt(rho_GND)`, with "
        f"`alpha = {ALPHA}`, `M = {TAYLOR_FACTOR}`, `G = {SHEAR_MODULUS_PA / 1e9:.1f} GPa`, "
        f"and `b = {BURGERS_VECTOR_M * 1e9:.3f} nm`.",
        "",
        "## Fit-free literature-parameter result",
        "",
        f"- Six-state MAPE: `{literature_metrics['mape_all_percent']:.2f}%`.",
        f"- Deformed-state MAPE: `{literature_metrics['mape_deformed_only_percent']:.2f}%`.",
        f"- Maximum absolute relative deviation: `{literature_metrics['max_absolute_relative_deviation_percent']:.2f}%`.",
        f"- RMSE: `{literature_metrics['rmse_MPa']:.1f} MPa`.",
        "",
        "## Bounded common-parameter calibration",
        "",
        "To correct the two identified transfer mismatches without state-specific fitting, "
        "the calibrated version replaces the reference-alloy matrix term by one common "
        "effective matrix stress and converts KAM-GND to an effective total dislocation density "
        "through one common density ratio. Solid solution and HAGB parameters remain fixed.",
        "",
        f"- Calibrated matrix term: `{calibration['sigma_matrix_calibrated_MPa']:.2f} MPa`.",
        f"- Dislocation strength multiplier: `{calibration['dislocation_strength_multiplier']:.3f}`.",
        f"- Effective `rho_total/rho_GND`: `{calibration['rho_total_to_gnd_ratio']:.3f}`.",
        f"- Six-state MAPE: `{calibrated_metrics['mape_all_percent']:.2f}%`.",
        f"- Deformed-state MAPE: `{calibrated_metrics['mape_deformed_only_percent']:.2f}%`.",
        f"- Maximum absolute relative deviation: `{calibrated_metrics['max_absolute_relative_deviation_percent']:.2f}%`.",
        f"- RMSE: `{calibrated_metrics['rmse_MPa']:.1f} MPa`.",
        "",
        "## Interpretation limits",
        "",
        "- The solid-solution term is constant because all states come from the same material batch. "
        "It explains the alloy-composition baseline but not the increase with cold reduction.",
        "- The O-equivalent literature model does not explicitly isolate the Fe contribution; therefore "
        "the numerical solid-solution term should be described as a literature-parameter estimate.",
        "- The project O-equivalent value (0.3484 wt.%) is slightly above the source model's reported "
        "0.14--0.32 wt.% validation interval. Its use is therefore a modest extrapolation, not an "
        "independently validated measurement of solid-solution strengthening.",
        "- KAM supplies an EBSD-derived GND estimate rather than total dislocation density. The model "
        "does not add a separate LAGB term, avoiding direct double counting of boundary-arranged dislocations.",
        "- The HAGB number-mean ECD changes little among the six states, so the Hall--Petch term is nearly constant.",
        "- No parameter is fitted to the six experimental yield strengths. The remaining underprediction "
        "after deformation may include statistically stored dislocations, residual stress, texture-dependent "
        "Taylor response, and uncertainty in effective strengthening length scale; these are not assigned quantitatively here.",
        "- The preceding statement applies to the literature-parameter version. The calibrated version "
        "uses two shared parameters fitted against all six yield strengths and is a descriptive closure, "
        "not an independent prediction or validation.",
        "",
        "## Sources",
        "",
        "- Solid-solution model: https://doi.org/10.1007/s12540-014-6004-8",
        "- Independent CP-Ti Hall--Petch parameter: https://doi.org/10.2320/matertrans.M2018033",
        "- Reference decomposition and matrix/Taylor parameters: attached Xia et al. reference paper.",
        "",
        "## Figure outputs",
        "",
        "- Six state figures: `figures/individual/` (editable SVG and 600 dpi PNG/TIFF).",
        "- Six-state montage: `figures/reference_three_factor_strengthening_montage` "
        "(editable SVG, 600 dpi PNG/TIFF, and vector PDF).",
        "- Calibrated six-state figures: `figures/calibrated/`, including individual figures "
        "and `reference_three_factor_strengthening_calibrated_montage` in SVG/PNG/TIFF/PDF.",
        "",
    ]
    (OUTPUT / "README.md").write_text("\n".join(report), encoding="utf-8")


def main() -> None:
    write_outputs(calculate())


if __name__ == "__main__":
    main()
