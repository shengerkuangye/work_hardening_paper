"""Create separate six-diameter OIM montages for IPF, PF, and ODF."""

from __future__ import annotations

import argparse
from pathlib import Path
from string import ascii_lowercase

from PIL import Image, ImageChops, ImageDraw

from assemble_oim_ipf_pf_odf_matrix import (
    SAMPLES,
    centered_lines,
    centered_text,
    draw_intensity_legend,
    font,
    parse_legend,
    prepare_panel,
)


PANEL_SIZE = 1024
OUTER_MARGIN = 50
COLUMN_GAP = 34
ROW_GAP = 34
PANEL_HEADER = 64


def panel_title(sample) -> str:
    return (
        f"{sample.diameter_mm:.2f} mm  |  "
        f"area reduction {sample.area_reduction_percent:.2f}%"
    )


def paste_and_verify(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    panel: Image.Image,
    x: int,
    y: int,
) -> None:
    canvas.paste(panel, (x, y))
    draw.rectangle(
        (x - 1, y - 1, x + PANEL_SIZE, y + PANEL_SIZE),
        outline="#666666",
        width=1,
    )
    pasted = canvas.crop((x, y, x + PANEL_SIZE, y + PANEL_SIZE))
    if ImageChops.difference(panel, pasted).getbbox() is not None:
        raise RuntimeError(f"Panel verification failed at ({x}, {y})")


def save_outputs(canvas: Image.Image, output_dir: Path, stem: str) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / f"{stem}.png"
    tiff_path = output_dir / f"{stem}.tif"
    canvas.save(png_path, format="PNG", dpi=(600, 600), optimize=True)
    canvas.save(tiff_path, format="TIFF", dpi=(600, 600), compression="tiff_lzw")
    return png_path, tiff_path


def draw_panel_heading(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    sample,
    panel_index: int,
    title_font,
    letter_font,
) -> None:
    draw.text((x + 2, y + 6), f"({ascii_lowercase[panel_index]})", font=letter_font, fill="black")
    centered_text(draw, x + PANEL_SIZE // 2, y + 8, panel_title(sample), title_font)


def build_ipf(source_root: Path, output_dir: Path) -> tuple[Path, Path]:
    ipf_dir = source_root / "oim_ipf"
    header_height = 390
    footer_height = 62
    cell_width = PANEL_SIZE
    row_height = PANEL_HEADER + PANEL_SIZE
    width = 2 * OUTER_MARGIN + 3 * cell_width + 2 * COLUMN_GAP
    height = header_height + 2 * row_height + ROW_GAP + footer_height
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)

    main_font = font(54, bold=True)
    sub_font = font(30)
    title_font = font(36)
    letter_font = font(39, bold=True)
    footer_font = font(26)

    centered_text(draw, width // 2, 12, "OIM Ti-Hex IPF-[001] maps", main_font)
    centered_lines(
        draw,
        width // 2,
        78,
        [
            "Color reference: sample A3 (map normal)",
            "Map horizontal direction: A2 = AD (bar axis)",
        ],
        sub_font,
        line_gap=4,
    )

    key_source = Image.open(ipf_dir / "tuli.jpg").convert("RGB")
    if key_source.size != (417, 433):
        key_source.close()
        raise ValueError(f"Unexpected IPF legend size: {key_source.size}")
    key_crop = key_source.crop((0, 228, 235, 386))
    key_source.close()
    key_display = key_crop.resize((330, 222), Image.Resampling.LANCZOS)
    key_crop.close()
    canvas.paste(key_display, (width // 2 - key_display.width // 2, 160))
    key_display.close()
    draw.line((OUTER_MARGIN, header_height - 1, width - OUTER_MARGIN, header_height - 1), fill="#888888", width=2)

    for index, sample in enumerate(SAMPLES):
        row, column = divmod(index, 3)
        cell_x = OUTER_MARGIN + column * (cell_width + COLUMN_GAP)
        row_y = header_height + row * (row_height + ROW_GAP)
        image_y = row_y + PANEL_HEADER
        draw_panel_heading(draw, cell_x, row_y, sample, index, title_font, letter_font)
        panel = prepare_panel(
            ipf_dir / f"oim_ipf_{sample.stem}.jpg",
            PANEL_SIZE,
            expect_size=(1024, 1024),
        )
        paste_and_verify(canvas, draw, panel, cell_x, image_y)
        panel.close()

    centered_text(
        draw,
        width // 2,
        height - 45,
        "IPF color identifies the Ti-Hex crystal direction parallel to sample A3; it is not an intensity scale.",
        footer_font,
        fill="#333333",
    )
    paths = save_outputs(canvas, output_dir, "oim_ipf_six_diameters")
    canvas.close()
    return paths


def build_pf(source_root: Path, output_dir: Path) -> tuple[Path, Path]:
    pf_dir = source_root / "oim_pf"
    header_height = 230
    footer_height = 72
    legend_width = 184
    cell_width = PANEL_SIZE + legend_width
    row_height = PANEL_HEADER + PANEL_SIZE
    width = 2 * OUTER_MARGIN + 3 * cell_width + 2 * COLUMN_GAP
    height = header_height + 2 * row_height + ROW_GAP + footer_height
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)

    main_font = font(54, bold=True)
    sub_font = font(30)
    title_font = font(36)
    letter_font = font(39, bold=True)
    legend_font = font(22)
    footer_font = font(26)

    centered_text(draw, width // 2, 12, "OIM Ti-Hex {0001} pole figures", main_font)
    centered_lines(
        draw,
        width // 2,
        80,
        [
            "c-axis distribution | top: +A1 | left: +A2 = AD | center: A3 (map normal)",
            "Harmonic series expansion: L=16, smoothing=5 deg, sample symmetry=Triclinic",
        ],
        sub_font,
        line_gap=8,
    )
    draw.line((OUTER_MARGIN, header_height - 1, width - OUTER_MARGIN, header_height - 1), fill="#888888", width=2)

    for index, sample in enumerate(SAMPLES):
        row, column = divmod(index, 3)
        cell_x = OUTER_MARGIN + column * (cell_width + COLUMN_GAP)
        row_y = header_height + row * (row_height + ROW_GAP)
        image_y = row_y + PANEL_HEADER
        draw_panel_heading(draw, cell_x, row_y, sample, index, title_font, letter_font)
        panel = prepare_panel(
            pf_dir / sample.pf_filename,
            PANEL_SIZE,
            expect_size=(1280, 1280),
        )
        paste_and_verify(canvas, draw, panel, cell_x, image_y)
        panel.close()
        bins = parse_legend(pf_dir / f"{sample.stem}.txt")
        draw_intensity_legend(
            draw,
            cell_x + PANEL_SIZE + 12,
            image_y + 112,
            bins,
            legend_font,
        )

    centered_text(
        draw,
        width // 2,
        height - 49,
        "Each panel retains its specimen-specific OIM auto bins; identical colors across panels are not identical intensity values.",
        footer_font,
        fill="#333333",
    )
    paths = save_outputs(canvas, output_dir, "oim_pf_0001_six_diameters")
    canvas.close()
    return paths


def build_odf(source_root: Path, output_dir: Path) -> tuple[Path, Path]:
    odf_dir = source_root / "oim_odf"
    header_height = 230
    footer_height = 72
    legend_width = 184
    cell_width = PANEL_SIZE + legend_width
    row_height = PANEL_HEADER + PANEL_SIZE
    width = 2 * OUTER_MARGIN + 3 * cell_width + 2 * COLUMN_GAP
    height = header_height + 2 * row_height + ROW_GAP + footer_height
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)

    main_font = font(54, bold=True)
    sub_font = font(30)
    title_font = font(36)
    letter_font = font(39, bold=True)
    legend_font = font(22)
    footer_font = font(26)

    centered_text(draw, width // 2, 12, "OIM Ti-Hex orientation distribution function sections", main_font)
    centered_lines(
        draw,
        width // 2,
        80,
        [
            "Bunge Euler space: phi1=0-90 deg, PHI=0-90 deg, phi2=0-60 deg (interval 10 deg)",
            "Harmonic series expansion: L=16, smoothing=5 deg, sample symmetry=Triclinic",
        ],
        sub_font,
        line_gap=8,
    )
    draw.line((OUTER_MARGIN, header_height - 1, width - OUTER_MARGIN, header_height - 1), fill="#888888", width=2)

    for index, sample in enumerate(SAMPLES):
        row, column = divmod(index, 3)
        cell_x = OUTER_MARGIN + column * (cell_width + COLUMN_GAP)
        row_y = header_height + row * (row_height + ROW_GAP)
        image_y = row_y + PANEL_HEADER
        draw_panel_heading(draw, cell_x, row_y, sample, index, title_font, letter_font)
        panel = prepare_panel(
            odf_dir / f"oim_odf_{sample.stem}.jpg",
            PANEL_SIZE,
            expect_size=(1024, 1024),
        )
        paste_and_verify(canvas, draw, panel, cell_x, image_y)
        panel.close()
        bins = parse_legend(odf_dir / f"{sample.stem}.txt")
        draw_intensity_legend(
            draw,
            cell_x + PANEL_SIZE + 12,
            image_y + 112,
            bins,
            legend_font,
        )

    centered_text(
        draw,
        width // 2,
        height - 49,
        "Each panel retains its specimen-specific OIM auto bins; compare peak locations directly, but compare intensities using the adjacent bins.",
        footer_font,
        fill="#333333",
    )
    paths = save_outputs(canvas, output_dir, "oim_odf_six_diameters")
    canvas.close()
    return paths


def parse_args() -> argparse.Namespace:
    repository = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-root",
        type=Path,
        default=repository / "data" / "ebsd_kpl_250221_7_df" / "images_oim",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repository / "results" / "oim_separate_montages",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    source = args.source_root.resolve()
    output = args.output_dir.resolve()
    results = {
        "IPF": build_ipf(source, output),
        "PF": build_pf(source, output),
        "ODF": build_odf(source, output),
    }
    for name, (png, tiff) in results.items():
        print(f"{name} PNG: {png}")
        print(f"{name} TIFF: {tiff}")
    print("Verified: each montage uses 7.00, 6.48, 6.02, 5.60, 5.25, 5.00 mm in row-major order")
