"""Extract OIM ODF sections and arrange them as diameter x phi2 (6 x 7)."""

from __future__ import annotations

import argparse
from pathlib import Path
from string import ascii_lowercase

from PIL import Image, ImageChops, ImageDraw

from assemble_oim_ipf_pf_odf_matrix import (
    SAMPLES,
    centered_lines,
    centered_text,
    font,
    interval_label,
    parse_legend,
)


PHI2_DEGREES = (0, 10, 20, 30, 40, 50, 60)
EXPECTED_SECTION_BOXES = (
    (27, 14, 315, 301),
    (369, 14, 656, 301),
    (710, 14, 998, 301),
    (27, 355, 315, 643),
    (369, 355, 656, 643),
    (710, 355, 998, 643),
    (27, 696, 315, 984),
)


def contiguous_runs(indices: list[int]) -> list[tuple[int, int]]:
    if not indices:
        return []
    runs: list[tuple[int, int]] = []
    start = previous = indices[0]
    for value in indices[1:]:
        if value != previous + 1:
            runs.append((start, previous + 1))
            start = value
        previous = value
    runs.append((start, previous + 1))
    return runs


def detect_section_boxes(image: Image.Image) -> list[tuple[int, int, int, int]]:
    """Detect the seven saturated square plots while excluding black angle text."""
    if image.size != (1024, 1024):
        raise ValueError(f"Expected 1024 x 1024 ODF source; found {image.size}")
    pixels = image.load()
    column_counts = [0] * image.width
    row_counts = [0] * image.height
    for y in range(image.height):
        for x in range(image.width):
            red, green, blue = pixels[x, y]
            if max(red, green, blue) - min(red, green, blue) > 25:
                column_counts[x] += 1
                row_counts[y] += 1

    x_runs = contiguous_runs([index for index, count in enumerate(column_counts) if count > 100])
    y_runs = contiguous_runs([index for index, count in enumerate(row_counts) if count > 100])
    if len(x_runs) != 3 or len(y_runs) != 3:
        raise ValueError(f"Expected 3 x-runs and 3 y-runs; found x={x_runs}, y={y_runs}")

    widths = [end - start for start, end in x_runs]
    heights = [end - start for start, end in y_runs]
    if max(widths) - min(widths) > 2 or max(heights) - min(heights) > 2:
        raise ValueError(f"Inconsistent section dimensions: widths={widths}, heights={heights}")
    if abs(round(sum(widths) / 3) - round(sum(heights) / 3)) > 1:
        raise ValueError(f"Detected ODF sections are not square: widths={widths}, heights={heights}")

    detected_boxes = [
        (x_runs[column][0], y_runs[row][0], x_runs[column][1], y_runs[row][1])
        for row, column in ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0))
    ]
    for detected, expected in zip(detected_boxes, EXPECTED_SECTION_BOXES, strict=True):
        if any(abs(left - right) > 2 for left, right in zip(detected, expected, strict=True)):
            raise ValueError(f"Detected section {detected} is inconsistent with expected box {expected}")
    return list(EXPECTED_SECTION_BOXES)


def draw_row_legend(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    row_height: int,
    bins,
    title_font,
    label_font,
) -> None:
    centered_lines(
        draw,
        x + 135,
        y + 5,
        ["OIM intensity", "1 = random | auto bins"],
        title_font,
        line_gap=2,
    )
    bar_x = x + 8
    bar_y = y + 68
    bar_width = 42
    block_height = (row_height - 82) // 8
    for display_row, bin_index in enumerate(reversed(range(8))):
        entry = bins[bin_index]
        y0 = bar_y + display_row * block_height
        y1 = y0 + block_height
        draw.rectangle((bar_x, y0, bar_x + bar_width, y1), fill=entry.color, outline="black", width=1)
        label = interval_label(bins, bin_index)
        box = draw.textbbox((0, 0), label, font=label_font)
        label_y = y0 + (block_height - (box[3] - box[1])) / 2 - box[1]
        draw.text((bar_x + bar_width + 10, label_y), label, font=label_font, fill="black")


def build_matrix(source_root: Path, output_dir: Path) -> tuple[Path, Path]:
    odf_dir = source_root / "oim_odf"
    cell_size = 420
    column_gap = 14
    row_gap = 18
    outer_margin = 48
    row_label_width = 285
    legend_width = 305
    header_height = 225
    footer_height = 70
    width = (
        2 * outer_margin
        + row_label_width
        + 7 * cell_size
        + 6 * column_gap
        + legend_width
    )
    height = header_height + 6 * cell_size + 5 * row_gap + footer_height
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)

    main_font = font(58, bold=True)
    subtitle_font = font(31)
    column_font = font(42, bold=True)
    row_font = font(43, bold=True)
    row_subfont = font(29)
    row_letter_font = font(38, bold=True)
    legend_title_font = font(27)
    legend_label_font = font(25)
    footer_font = font(28)

    grid_x = outer_margin + row_label_width
    legend_x = grid_x + 7 * cell_size + 6 * column_gap
    centered_text(draw, width // 2, 10, "OIM Ti-Hex ODF sections: diameter × φ2", main_font)
    centered_lines(
        draw,
        width // 2,
        79,
        [
            "Bunge Euler space | horizontal: φ1 = 0–90° | vertical: Φ = 0–90°",
            "Harmonic series expansion: L=16, smoothing=5°, sample symmetry=Triclinic",
        ],
        subtitle_font,
        line_gap=5,
    )
    centered_text(draw, outer_margin + row_label_width // 2, 174, "Diameter", column_font)
    for column, phi2 in enumerate(PHI2_DEGREES):
        center_x = grid_x + column * (cell_size + column_gap) + cell_size // 2
        centered_text(draw, center_x, 174, f"φ2 = {phi2}°", column_font)
    centered_text(draw, legend_x + legend_width // 2, 174, "Intensity bins", column_font)
    draw.line((outer_margin, header_height - 1, width - outer_margin, header_height - 1), fill="#777777", width=2)

    reference_boxes: list[tuple[int, int, int, int]] | None = None
    pasted_sections: list[tuple[Image.Image, int, int]] = []
    for row, sample in enumerate(SAMPLES):
        row_y = header_height + row * (cell_size + row_gap)
        source_path = odf_dir / f"oim_odf_{sample.stem}.jpg"
        source = Image.open(source_path).convert("RGB")
        boxes = detect_section_boxes(source)
        if reference_boxes is None:
            reference_boxes = boxes
        elif boxes != reference_boxes:
            source.close()
            raise ValueError(f"ODF section boxes differ for {source_path}: {boxes} vs {reference_boxes}")

        draw.text((outer_margin + 2, row_y + 8), f"({ascii_lowercase[row]})", font=row_letter_font, fill="black")
        centered_text(
            draw,
            outer_margin + row_label_width // 2,
            row_y + cell_size // 2 - 65,
            f"{sample.diameter_mm:.2f} mm",
            row_font,
        )
        centered_text(
            draw,
            outer_margin + row_label_width // 2,
            row_y + cell_size // 2,
            "area reduction",
            row_subfont,
            fill="#333333",
        )
        centered_text(
            draw,
            outer_margin + row_label_width // 2,
            row_y + cell_size // 2 + 42,
            f"{sample.area_reduction_percent:.2f}%",
            row_subfont,
            fill="#333333",
        )

        for column, box in enumerate(boxes):
            section = source.crop(box).resize((cell_size, cell_size), Image.Resampling.LANCZOS)
            x = grid_x + column * (cell_size + column_gap)
            canvas.paste(section, (x, row_y))
            draw.rectangle((x - 1, row_y - 1, x + cell_size, row_y + cell_size), outline="#555555", width=1)
            pasted_sections.append((section, x, row_y))
        source.close()

        bins = parse_legend(odf_dir / f"{sample.stem}.txt")
        draw_row_legend(
            draw,
            legend_x + 8,
            row_y,
            cell_size,
            bins,
            legend_title_font,
            legend_label_font,
        )
        if row < len(SAMPLES) - 1:
            separator_y = row_y + cell_size + row_gap // 2
            draw.line((outer_margin, separator_y, width - outer_margin, separator_y), fill="#D0D0D0", width=1)

    for section, x, y in pasted_sections:
        pasted = canvas.crop((x, y, x + cell_size, y + cell_size))
        if ImageChops.difference(section, pasted).getbbox() is not None:
            raise RuntimeError(f"ODF section verification failed at ({x}, {y})")
        section.close()

    centered_text(
        draw,
        width // 2,
        height - 47,
        "Each row retains its specimen-specific OIM auto bins; compare locations directly and use the row legend for intensity.",
        footer_font,
        fill="#333333",
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / "oim_odf_6x7_diameter_phi2_matrix.png"
    tiff_path = output_dir / "oim_odf_6x7_diameter_phi2_matrix.tif"
    canvas.save(png_path, format="PNG", dpi=(600, 600), optimize=True)
    canvas.save(tiff_path, format="TIFF", dpi=(600, 600), compression="tiff_lzw")

    # Keep the previously delivered ODF filenames current, so they no longer
    # point to the superseded 2 x 3 whole-image montage.
    canvas.save(output_dir / "oim_odf_six_diameters.png", format="PNG", dpi=(600, 600), optimize=True)
    canvas.save(
        output_dir / "oim_odf_six_diameters.tif",
        format="TIFF",
        dpi=(600, 600),
        compression="tiff_lzw",
    )
    canvas.close()

    print(f"Detected source crop boxes: {reference_boxes}")
    return png_path, tiff_path


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
    png, tiff = build_matrix(args.source_root.resolve(), args.output_dir.resolve())
    print(f"Created PNG: {png}")
    print(f"Created TIFF: {tiff}")
    print("Rows: 7.00, 6.48, 6.02, 5.60, 5.25, 5.00 mm")
    print("Columns: phi2 = 0, 10, 20, 30, 40, 50, 60 deg")
    print("Verified: 42 extracted ODF sections, no rotation, clean degree labels in the matrix header")
