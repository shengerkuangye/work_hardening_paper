#!/usr/bin/env python3
"""Calculate literature-form and closest-fit LAGB strengthening models.

The script reads existing derived EBSD tables and the laboratory tensile
summary.  It never modifies raw experimental data.  Outputs are JSON and a
Markdown audit report under results/lagb_strengthening_models/.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TENSILE = ROOT / "data/tensile_data/gr4b23271_lab_tensile_by_diameter.csv"
HAGB = ROOT / "results/mtex_grain_size_distribution_matrix/hagb_grain_size_summary.csv"
GB_SUMMARY = ROOT / "results/mtex_grain_boundary_misorientation/grain_boundary_misorientation_summary.csv"
GB_DISTRIBUTION = ROOT / "results/mtex_grain_boundary_misorientation/grain_boundary_misorientation_distribution.csv"
OUTPUT_DIR = ROOT / "results/lagb_strengthening_models"

SAMPLES = ["7d", "6.48d", "6.02d", "5.6d", "5.25d", "5d"]
TENSILE_DIAMETERS = [7.0, 6.5, 6.02, 5.6, 5.25, 5.0]

# Independent parameters for commercially pure alpha-Ti.
K_HAGB_MPA_UM_HALF = 330.0
TAYLOR_FACTOR = 3.0
ALPHA = 0.2
SHEAR_MODULUS_MPA = 45_600.0
BURGERS_VECTOR_UM = 0.000295


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def keyed(rows: list[dict[str, str]], field: str) -> dict[str, dict[str, str]]:
    return {row[field]: row for row in rows}


def weighted_median(values: list[float], weights: list[float]) -> float:
    ordered = sorted(zip(values, weights))
    half = sum(weights) / 2.0
    cumulative = 0.0
    for value, weight in ordered:
        cumulative += weight
        if cumulative >= half:
            return value
    return ordered[-1][0]


def metrics(rows: list[dict[str, float]], prediction_key: str) -> dict[str, float]:
    errors = [
        100.0 * (row[prediction_key] - row["yield_experimental_MPa"])
        / row["yield_experimental_MPa"]
        for row in rows
    ]
    residuals = [
        row[prediction_key] - row["yield_experimental_MPa"] for row in rows
    ]
    return {
        "mape_all_percent": sum(abs(value) for value in errors) / len(errors),
        "mape_deformed_only_percent": sum(abs(value) for value in errors[1:])
        / (len(errors) - 1),
        "max_absolute_relative_deviation_percent": max(abs(value) for value in errors),
        "rmse_MPa": math.sqrt(sum(value * value for value in residuals) / len(residuals)),
    }


def calculate() -> dict[str, object]:
    tensile_by_diameter = {
        float(row["nominal_diameter_group_mm"]): row for row in read_rows(TENSILE)
    }
    hagb_by_sample = keyed(read_rows(HAGB), "sample")
    gb_by_sample = keyed(read_rows(GB_SUMMARY), "sample")

    theta_numerator = {sample: 0.0 for sample in SAMPLES}
    theta_denominator = {sample: 0.0 for sample in SAMPLES}
    for row in read_rows(GB_DISTRIBUTION):
        sample = row["sample"]
        if sample not in theta_numerator:
            continue
        center_deg = float(row["bin_center_deg"])
        if 1.0 <= center_deg < 15.0:
            length_um = float(row["boundary_length_um"])
            theta_numerator[sample] += center_deg * length_um
            theta_denominator[sample] += length_um

    rows: list[dict[str, float | str]] = []
    for sample, tensile_diameter in zip(SAMPLES, TENSILE_DIAMETERS):
        tensile = tensile_by_diameter[tensile_diameter]
        hagb = hagb_by_sample[sample]
        gb = gb_by_sample[sample]
        theta_deg = theta_numerator[sample] / theta_denominator[sample]
        theta_rad = math.radians(theta_deg)
        lagb_length_density = float(gb["lagb_length_density_um_per_um2"])

        # For an isotropic two-dimensional line network, intersections per
        # unit test-line length are 2*lambda/pi, hence L = pi/(2*lambda).
        lagb_intercept_um = math.pi / (2.0 * lagb_length_density)
        hagb_size_um = float(hagb["ecd_number_mean_um"])
        sigma_hagb = K_HAGB_MPA_UM_HALF / math.sqrt(hagb_size_um)
        sigma_lagb = (
            TAYLOR_FACTOR
            * ALPHA
            * SHEAR_MODULUS_MPA
            * math.sqrt(
                3.0 * BURGERS_VECTOR_UM * theta_rad / lagb_intercept_um
            )
        )
        rows.append(
            {
                "sample": sample,
                # Preserve the EBSD state label (6.48 mm/14.31%) while mapping
                # it to the laboratory tensile group reported as nominal 6.5 mm.
                "diameter_mm": float(hagb["diameter_mm"]),
                "cold_reduction_percent": float(hagb["cold_reduction_percent"]),
                "yield_experimental_MPa": float(tensile["Rp0.2_MPa_mean"]),
                "yield_std_MPa": float(tensile["Rp0.2_MPa_std"]),
                "hagb_number_mean_ecd_um": hagb_size_um,
                "lagb_length_density_um_per_um2": lagb_length_density,
                "lagb_length_weighted_mean_misorientation_deg": theta_deg,
                "lagb_stereological_intercept_um": lagb_intercept_um,
                "sigma_hagb_MPa": sigma_hagb,
                "sigma_lagb_literature_MPa": sigma_lagb,
            }
        )

    baseline = (
        rows[0]["yield_experimental_MPa"]
        - rows[0]["sigma_hagb_MPa"]
        - rows[0]["sigma_lagb_literature_MPa"]
    )
    for row in rows:
        row["yield_strict_MPa"] = (
            baseline + row["sigma_hagb_MPa"] + row["sigma_lagb_literature_MPa"]
        )
        row["strict_relative_deviation_percent"] = 100.0 * (
            row["yield_strict_MPa"] - row["yield_experimental_MPa"]
        ) / row["yield_experimental_MPa"]

    first = rows[0]
    q_breakpoints: list[float] = []
    q_weights: list[float] = []
    for row in rows[1:]:
        base_curve = (
            first["yield_experimental_MPa"]
            + row["sigma_hagb_MPa"]
            - first["sigma_hagb_MPa"]
        )
        x = row["sigma_lagb_literature_MPa"] - first["sigma_lagb_literature_MPa"]
        q_breakpoints.append((row["yield_experimental_MPa"] - base_curve) / x)
        q_weights.append(abs(x) / row["yield_experimental_MPa"])

    # This weighted median is the exact minimizer of mean absolute relative
    # error for the one-parameter model.
    q_closest = weighted_median(q_breakpoints, q_weights)
    closest_baseline = (
        first["yield_experimental_MPa"]
        - first["sigma_hagb_MPa"]
        - q_closest * first["sigma_lagb_literature_MPa"]
    )
    for row in rows:
        row["sigma_lagb_closest_MPa"] = q_closest * row["sigma_lagb_literature_MPa"]
        row["yield_closest_MPa"] = (
            closest_baseline
            + row["sigma_hagb_MPa"]
            + row["sigma_lagb_closest_MPa"]
        )
        row["closest_relative_deviation_percent"] = 100.0 * (
            row["yield_closest_MPa"] - row["yield_experimental_MPa"]
        ) / row["yield_experimental_MPa"]

    loo_errors: list[float] = []
    for omitted in range(1, len(rows)):
        breakpoints = [
            q_breakpoints[index - 1]
            for index in range(1, len(rows))
            if index != omitted
        ]
        weights = [
            q_weights[index - 1]
            for index in range(1, len(rows))
            if index != omitted
        ]
        q_loo = weighted_median(breakpoints, weights)
        row = rows[omitted]
        prediction = (
            first["yield_experimental_MPa"]
            + row["sigma_hagb_MPa"]
            - first["sigma_hagb_MPa"]
            + q_loo
            * (
                row["sigma_lagb_literature_MPa"]
                - first["sigma_lagb_literature_MPa"]
            )
        )
        loo_errors.append(
            100.0 * (prediction - row["yield_experimental_MPa"])
            / row["yield_experimental_MPa"]
        )

    strict_metrics = metrics(rows, "yield_strict_MPa")
    closest_metrics = metrics(rows, "yield_closest_MPa")
    closest_metrics["leave_one_deformed_state_out_mape_percent"] = sum(
        abs(value) for value in loo_errors
    ) / len(loo_errors)
    closest_metrics["leave_one_deformed_state_out_max_percent"] = max(
        abs(value) for value in loo_errors
    )

    return {
        "model_scope": (
            "HAGB Hall-Petch plus LAGB dislocation-boundary strengthening; "
            "KAM-GND is not added separately to avoid direct double counting"
        ),
        "source_urls": [
            "https://doi.org/10.2320/matertrans.M2018033",
            "https://doi.org/10.1016/S1359-6454(97)00365-0",
        ],
        "parameters": {
            "k_hagb_MPa_um_half": K_HAGB_MPA_UM_HALF,
            "M": TAYLOR_FACTOR,
            "alpha": ALPHA,
            "G_MPa": SHEAR_MODULUS_MPA,
            "b_um": BURGERS_VECTOR_UM,
            "strict_baseline_MPa": baseline,
            "closest_q": q_closest,
            "closest_baseline_MPa": closest_baseline,
            "closest_equivalent_alpha_if_M_fixed": ALPHA * q_closest,
            "closest_equivalent_M_if_alpha_fixed": TAYLOR_FACTOR * q_closest,
            "closest_equivalent_intercept_ratio_if_only_L_changes": 1.0
            / (q_closest * q_closest),
        },
        "strict_metrics": strict_metrics,
        "closest_metrics": closest_metrics,
        "rows": rows,
    }


def write_report(result: dict[str, object]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT_DIR / "lagb_strengthening_models.json"
    json_path.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    params = result["parameters"]
    strict = result["strict_metrics"]
    closest = result["closest_metrics"]
    lines = [
        "# LAGB strengthening models",
        "",
        "## Model scope",
        "",
        "The calculation follows the literature decomposition "
        "`sigma_y = sigma_0 + sigma_HAGB + sigma_LAGB`. KAM-derived GND "
        "strengthening is not added separately because the LAGB model already "
        "represents dislocations arranged in low-angle boundaries.",
        "",
        "Strict model: literature-form equation and independent CP-alpha-Ti "
        "parameters; the 7 mm condition calibrates the matrix/base term only.",
        "",
        "Closest model: the same structure, with one dimensionless multiplier "
        "`q` on the LAGB term chosen to minimize MAPE over the five deformed "
        "conditions.",
        "",
        "## Parameters and metrics",
        "",
        f"- Strict baseline: {params['strict_baseline_MPa']:.3f} MPa",
        f"- Closest q: {params['closest_q']:.6f}",
        f"- Closest baseline: {params['closest_baseline_MPa']:.3f} MPa",
        f"- Strict MAPE, deformed only: {strict['mape_deformed_only_percent']:.3f}%",
        f"- Closest MAPE, deformed only: {closest['mape_deformed_only_percent']:.3f}%",
        f"- Closest leave-one-deformed-state-out MAPE: "
        f"{closest['leave_one_deformed_state_out_mape_percent']:.3f}%",
        "",
        "## Results",
        "",
        "| Diameter (mm) | Experimental (MPa) | HAGB (MPa) | LAGB strict (MPa) | Strict prediction (MPa) | Strict deviation (%) | Closest prediction (MPa) | Closest deviation (%) |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in result["rows"]:
        lines.append(
            f"| {row['diameter_mm']:.2f} | {row['yield_experimental_MPa']:.1f} | "
            f"{row['sigma_hagb_MPa']:.1f} | "
            f"{row['sigma_lagb_literature_MPa']:.1f} | "
            f"{row['yield_strict_MPa']:.1f} | "
            f"{row['strict_relative_deviation_percent']:+.2f} | "
            f"{row['yield_closest_MPa']:.1f} | "
            f"{row['closest_relative_deviation_percent']:+.2f} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "The strict version is the defensible literature-form estimate with "
            "currently available EBSD quantities. Its stereological intercept "
            "assumes a two-dimensional isotropic boundary network: "
            "`L = pi/(2 lambda_LAGB)`. Direct directional line-intercept "
            "measurements should replace this approximation for the elongated "
            "cold-worked structure.",
            "",
            "The closest version is a fitted effective model, not independent "
            "validation. The multiplier can represent boundary anisotropy, "
            "line-intercept conversion, alpha, Taylor-factor, and unresolved "
            "substructure effects; it must not be assigned uniquely to any one "
            "mechanism without additional measurements.",
            "",
            "Sources:",
            "",
            "- https://doi.org/10.2320/matertrans.M2018033",
            "- https://doi.org/10.1016/S1359-6454(97)00365-0",
            "",
        ]
    )
    (OUTPUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    write_report(calculate())
