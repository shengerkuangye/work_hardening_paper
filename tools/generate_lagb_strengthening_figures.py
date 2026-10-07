#!/usr/bin/env python3
"""Generate Fig. 15-style HAGB/LAGB strengthening figures.

Creates one figure per deformation state for the strict and closest models,
plus one six-panel montage for each model. SVG is the editable master; PNG,
TIFF and PDF are exported at 600 dpi or as vector output.
"""

from __future__ import annotations

import html
import json
import subprocess
from pathlib import Path

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results/lagb_strengthening_models/lagb_strengthening_models.json"
OUTPUT = ROOT / "results/lagb_strengthening_models/figures"
FONT = "Liberation Serif"
Y_MAX = 1100.0

COLORS = {
    "base": ("#f07d78", "#ffe4e1"),
    "hagb": ("#58d66d", "#dcffdf"),
    "lagb": ("#e8d95c", "#fff8c7"),
    "experimental": ("#3aa5e8", "#d4efff"),
}


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def defs(prefix: str) -> str:
    blocks = []
    for name, (strong, pale) in COLORS.items():
        blocks.append(
            f'<linearGradient id="{prefix}-{name}" x1="0" y1="0" x2="1" y2="0">'
            f'<stop offset="0" stop-color="{pale}"/>'
            f'<stop offset="0.5" stop-color="{strong}"/>'
            f'<stop offset="1" stop-color="{pale}"/>'
            "</linearGradient>"
        )
    return "".join(blocks)


def text(x: float, y: float, value: str, size: float, *, anchor: str = "middle", weight: str = "normal", fill: str = "#111111", rotate: float | None = None) -> str:
    transform = f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ""
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" '
        f'font-family="{FONT}" font-size="{size:.2f}" font-weight="{weight}" '
        f'fill="{fill}"{transform}>{esc(value)}</text>'
    )


def legend(prefix: str, x: float, y: float, width: float, size: float = 30.0) -> str:
    items = [
        ("base", "σbase"),
        ("hagb", "σHAGB"),
        ("lagb", "σLAGB"),
        ("experimental", "σexperimental"),
    ]
    item_width = width / len(items)
    parts = []
    for index, (key, label) in enumerate(items):
        left = x + index * item_width
        parts.append(
            f'<rect x="{left:.2f}" y="{y:.2f}" width="{size * 1.35:.2f}" '
            f'height="{size * 0.7:.2f}" fill="url(#{prefix}-{key})"/>'
        )
        parts.append(
            text(left + size * 1.55, y + size * 0.58, label, size * 0.72, anchor="start")
        )
    return "".join(parts)


def panel(
    row: dict[str, float | str],
    mode: str,
    prefix: str,
    *,
    width: float = 1200.0,
    height: float = 1100.0,
    include_legend: bool = True,
    panel_label: str | None = None,
) -> str:
    left, right = 150.0, width - 65.0
    top = 190.0 if include_legend else 115.0
    bottom = height - 155.0
    plot_height = bottom - top
    plot_width = right - left
    calc_x = left + plot_width * 0.36
    exp_x = left + plot_width * 0.72
    bar_width = plot_width * 0.19
    base = float(row["yield_strict_MPa"]) - float(row["sigma_hagb_MPa"]) - float(row["sigma_lagb_literature_MPa"])
    if mode == "closest":
        base = float(row["yield_closest_MPa"]) - float(row["sigma_hagb_MPa"]) - float(row["sigma_lagb_closest_MPa"])
        lagb = float(row["sigma_lagb_closest_MPa"])
        prediction = float(row["yield_closest_MPa"])
        deviation = float(row["closest_relative_deviation_percent"])
    else:
        lagb = float(row["sigma_lagb_literature_MPa"])
        prediction = float(row["yield_strict_MPa"])
        deviation = float(row["strict_relative_deviation_percent"])
    hagb = float(row["sigma_hagb_MPa"])
    experimental = float(row["yield_experimental_MPa"])

    def sy(value: float) -> float:
        return bottom - value / Y_MAX * plot_height

    parts = []
    if panel_label:
        parts.append(text(35, 45, panel_label, 38, anchor="start", weight="bold"))
    title = f"{float(row['diameter_mm']):.2f} mm   ({float(row['cold_reduction_percent']):.2f}%)"
    parts.append(text(width / 2, 58, title, 38, weight="bold"))
    if include_legend:
        parts.append(legend(prefix, left + 10, 94, plot_width - 20, 32))

    # Frame and ticks.
    parts.append(
        f'<rect x="{left:.2f}" y="{top:.2f}" width="{plot_width:.2f}" height="{plot_height:.2f}" '
        'fill="#ffffff" stroke="#111111" stroke-width="4"/>'
    )
    for tick in range(0, 1101, 200):
        y = sy(float(tick))
        parts.append(f'<line x1="{left - 12:.2f}" y1="{y:.2f}" x2="{left:.2f}" y2="{y:.2f}" stroke="#111" stroke-width="4"/>')
        parts.append(text(left - 22, y + 10, str(tick), 30, anchor="end"))
    parts.append(text(46, (top + bottom) / 2, "Yield strength (MPa)", 38, weight="bold", rotate=-90))

    # Stacked calculated bar.
    cumulative = 0.0
    segments = [("base", base), ("hagb", hagb), ("lagb", lagb)]
    for key, value in segments:
        y_top = sy(cumulative + value)
        y_bottom = sy(cumulative)
        parts.append(
            f'<rect x="{calc_x - bar_width / 2:.2f}" y="{y_top:.2f}" width="{bar_width:.2f}" '
            f'height="{y_bottom - y_top:.2f}" fill="url(#{prefix}-{key})"/>'
        )
        parts.append(text(calc_x, (y_top + y_bottom) / 2 + 9, f"{value:.1f}", 28, weight="bold"))
        cumulative += value
    parts.append(
        f'<rect x="{calc_x - bar_width / 2:.2f}" y="{sy(prediction):.2f}" width="{bar_width:.2f}" '
        f'height="{bottom - sy(prediction):.2f}" fill="none" stroke="#222" stroke-width="2"/>'
    )

    # Experimental bar.
    parts.append(
        f'<rect x="{exp_x - bar_width / 2:.2f}" y="{sy(experimental):.2f}" width="{bar_width:.2f}" '
        f'height="{bottom - sy(experimental):.2f}" fill="url(#{prefix}-experimental)" stroke="#228ac8" stroke-width="2"/>'
    )
    parts.append(text(exp_x, sy(experimental / 2) + 11, f"{experimental:.1f}", 32, weight="bold", fill="#ffffff"))

    # Values and category labels.
    parts.append(text(calc_x, sy(prediction) - 32, f"{prediction:.1f} MPa", 30, weight="bold"))
    parts.append(text(calc_x, sy(prediction) - 3, f"δ = {deviation:+.2f}%", 25))
    parts.append(text(calc_x, bottom + 50, "Calculated", 32, weight="bold"))
    parts.append(text(exp_x, bottom + 50, "Experimental", 32, weight="bold"))
    return "".join(parts)


def svg_document(body: str, width_in: float, height_in: float, view_width: int, view_height: int, prefix: str, title_value: str) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width_in}in" height="{height_in}in" '
        f'viewBox="0 0 {view_width} {view_height}">'
        f'<title>{esc(title_value)}</title><defs>{defs(prefix)}</defs>'
        '<rect width="100%" height="100%" fill="#ffffff"/>'
        f'{body}</svg>\n'
    )


def export(svg_path: Path, dpi: int = 600, *, include_pdf: bool = True) -> None:
    png_path = svg_path.with_suffix(".png")
    pdf_path = svg_path.with_suffix(".pdf")
    tif_path = svg_path.with_suffix(".tif")
    subprocess.run(
        ["inkscape", str(svg_path), "--export-type=png", f"--export-dpi={dpi}", f"--export-filename={png_path}"],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if include_pdf:
        subprocess.run(
            ["inkscape", str(svg_path), "--export-type=pdf", f"--export-filename={pdf_path}"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    with Image.open(png_path) as image:
        image.convert("RGB").save(tif_path, compression="tiff_lzw", dpi=(dpi, dpi))


def generate_individual(rows: list[dict[str, float | str]], mode: str) -> list[dict[str, str]]:
    directory = OUTPUT / mode / "individual"
    directory.mkdir(parents=True, exist_ok=True)
    manifest = []
    for row in rows:
        sample = str(row["sample"])
        prefix = f"{mode}-{sample}"
        body = panel(row, mode, prefix)
        svg = svg_document(body, 4.0, 4.0, 1200, 1100, prefix, f"{mode} strengthening - {sample}")
        svg_path = directory / f"{mode}_{sample}.svg"
        svg_path.write_text(svg, encoding="utf-8")
        export(svg_path, include_pdf=False)
        manifest.append({"sample": sample, "svg": str(svg_path.relative_to(ROOT))})
    return manifest


def generate_montage(rows: list[dict[str, float | str]], mode: str) -> Path:
    width, height = 3600, 2200
    prefix = f"{mode}-montage"
    title_value = "Literature-parameter strengthening model" if mode == "strict" else "One-parameter calibrated strengthening model"
    parts = [text(width / 2, 72, title_value, 50, weight="bold")]
    parts.append(legend(prefix, 710, 102, 2180, 38))
    labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
    cell_width, cell_height = 1200, 980
    for index, row in enumerate(rows):
        col, line = index % 3, index // 3
        local = panel(row, mode, prefix, width=1200, height=1100, include_legend=False, panel_label=labels[index])
        parts.append(
            f'<g transform="translate({col * cell_width},{180 + line * cell_height}) scale(1,0.89)">{local}</g>'
        )
    svg = svg_document("".join(parts), 12.0, 7.35, width, height, prefix, title_value)
    svg_path = OUTPUT / f"{mode}_strengthening_montage.svg"
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    svg_path.write_text(svg, encoding="utf-8")
    export(svg_path)
    return svg_path


def main() -> None:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    rows = data["rows"]
    manifest = {"source": str(INPUT.relative_to(ROOT)), "models": {}}
    for mode in ("strict", "closest"):
        individuals = generate_individual(rows, mode)
        montage = generate_montage(rows, mode)
        manifest["models"][mode] = {
            "individual": individuals,
            "montage": str(montage.relative_to(ROOT)),
        }
    (OUTPUT / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
