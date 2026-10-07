"""Assemble manuscript statistics from existing summaries; never change raw data.

Use curve-based strength after the recorded user exclusion, and laboratory
post-fracture A/Z records. The 6.48 legacy state is named 6.5 mm in the paper.
No curves, figures, yield fits, or strengthening models are regenerated.
"""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CURVE_DIR = ROOT / "data/tensile_data/gr4b23271_cold_deformation_2_raw_csv"
CURVE_SOURCE = CURVE_DIR / "gr4b23271_tensile_final_by_diameter_user_exclude.csv"
STATUS_SOURCE = CURVE_DIR / "gr4b23271_cold_deformation_tensile_summary.csv"
EXCLUSION_SOURCE = CURVE_DIR / "gr4b23271_tensile_user_exclude_detail.csv"
LAB_SOURCE = ROOT / "data/tensile_data/gr4b23271_lab_tensile_by_diameter.csv"
OUTPUT = ROOT / "tables/gr4b23271_mechanical_properties_curve_basis.csv"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def assemble() -> list[dict[str, object]]:
    lab = {float(r["nominal_diameter_group_mm"]): r for r in read_rows(LAB_SOURCE)}
    statuses = {r["sample"]: r["yield_status"] for r in read_rows(STATUS_SOURCE)}
    decisions = {r["sample"]: r for r in read_rows(EXCLUSION_SOURCE)}
    assembled = []
    curves = read_rows(CURVE_SOURCE)
    baseline = next(r for r in curves if float(r["diameter_mm"]) == 7)
    for curve in curves:
        legacy_diameter = float(curve["diameter_mm"])
        nominal_diameter = 6.5 if legacy_diameter == 6.48 else legacy_diameter
        laboratory = lab[nominal_diameter]
        included = [s for s in curve["included_samples"].split(";") if s]
        excluded = [s for s in curve["excluded_samples"].split(";") if s]
        n_curve = int(curve["valid_repeats_n"])
        if len(included) != n_curve:
            raise ValueError(f"Repeat-count mismatch: {legacy_diameter}")
        if any(statuses[s] != "ok" for s in included):
            raise ValueError(f"Unreliable yield value included: {legacy_diameter}")
        if any(decisions[s]["include_in_final_average"] != "yes" for s in included):
            raise ValueError(f"Inclusion differs from recorded decision: {legacy_diameter}")
        if any(decisions[s]["include_in_final_average"] != "no" for s in excluded):
            raise ValueError(f"Exclusion differs from recorded decision: {legacy_diameter}")
        reduction = 100.0 * (1.0 - (nominal_diameter / 7.0) ** 2)
        row = {
            "sample_name": (
                "M-7.00" if nominal_diameter == 7
                else "Y-6.5" if nominal_diameter == 6.5
                else f"Y-{nominal_diameter:.2f}"
            ),
            "nominal_diameter_mm": f"{nominal_diameter:g}",
            "nominal_area_reduction_percent": f"{reduction:.8f}",
            "legacy_curve_diameter_mm": curve["diameter_mm"],
            "legacy_curve_reduction_percent": curve["cold_reduction_percent_reference"],
            "lab_measured_diameter_mean_mm": laboratory["measured_d0_mm_mean"],
            "curve_n_valid": n_curve,
            "curve_n_total": curve["total_repeats_n"],
            "included_curve_samples": curve["included_samples"],
            "excluded_curve_samples": curve["excluded_samples"],
            "excluded_original_yield_status": ";".join(statuses[s] for s in excluded),
            "excluded_recorded_reason": ";".join(decisions[s]["exclude_reason"] for s in excluded),
            "Rp0.2_MPa_mean": curve["Rp0.2_MPa_mean"],
            "Rp0.2_MPa_sd": curve["Rp0.2_MPa_std"] if n_curve > 1 else "",
            "Rm_MPa_mean": curve["UTS_engineering_MPa_mean"],
            "Rm_MPa_sd": curve["UTS_engineering_MPa_std"] if n_curve > 1 else "",
            "engineering_strain_at_Rm_percent_mean": curve["uniform_elongation_percent_mean"],
            "engineering_strain_at_Rm_percent_sd": curve["uniform_elongation_percent_std"] if n_curve > 1 else "",
            "post_fracture_n_lab": laboratory["n"],
            "A_percent_mean": laboratory["total_elongation_A_percent_mean"],
            "A_percent_sd": laboratory["total_elongation_A_percent_std"],
            "Z_percent_mean": laboratory["reduction_of_area_Z_percent_mean"],
            "Z_percent_sd": laboratory["reduction_of_area_Z_percent_std"],
            "curve_statistics_source": str(CURVE_SOURCE.relative_to(ROOT)),
            "curve_status_source": str(STATUS_SOURCE.relative_to(ROOT)),
            "exclusion_source": str(EXCLUSION_SOURCE.relative_to(ROOT)),
            "post_fracture_source": str(LAB_SOURCE.relative_to(ROOT)),
            "strain_definition": "engineering strain at maximum engineering stress; not Ag or post-fracture A",
            "single_curve_sd_policy": "blank when n=1; no estimate of repeat scatter",
            "post_fracture_policy": "laboratory A/Z kept for both recorded repeats; curve exclusion is not applied to A/Z",
        }
        rp = float(curve["Rp0.2_MPa_mean"])
        rm = float(curve["UTS_engineering_MPa_mean"])
        row.update({
            "Rp0.2_change_vs_initial_percent": f"{100 * (rp / float(baseline['Rp0.2_MPa_mean']) - 1):.8f}",
            "Rm_change_vs_initial_percent": f"{100 * (rm / float(baseline['UTS_engineering_MPa_mean']) - 1):.8f}",
            "A_change_vs_initial_percent": f"{100 * (float(row['A_percent_mean']) / float(lab[7]['total_elongation_A_percent_mean']) - 1):.8f}",
            "Z_change_vs_initial_percent": f"{100 * (float(row['Z_percent_mean']) / float(lab[7]['reduction_of_area_Z_percent_mean']) - 1):.8f}",
            "Rp0.2_over_Rm_ratio_of_means": f"{rp / rm:.8f}",
            "Rm_minus_Rp0.2_difference_of_means_MPa": f"{rm - rp:.8f}",
        })
        assembled.append(row)
    if [r["nominal_diameter_mm"] for r in assembled] != ["7", "6.5", "6.02", "5.6", "5.25", "5"]:
        raise ValueError("Unexpected sample grouping or order")
    return assembled


if __name__ == "__main__":
    rows = assemble()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} states to {OUTPUT}")
