"""Create a 6 x 3 OIM IPF/PF/ODF comparison matrix.

Rows are ordered by decreasing specimen diameter.  Each PF and ODF panel is
paired with its own OIM-exported TXT intensity bins because OIM used an
independent automatic scale for every specimen.  Source files are read-only.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
from string import ascii_lowercase

from PIL import Image, ImageChops, ImageDraw, ImageFont


@dataclass(frozen=True)
class Sample:
    stem: str
    diameter_mm: float
    pf_filename: str

    @property
    def area_reduction_percent(self) -> float:
        return 100.0 * (1.0 - (self.diameter_mm / 7.0) ** 2)


@dataclass(frozen=True)
class LegendBin:
    operator: str
    threshold: float
    color: tuple[int, int, int]


SAMPLES = (
    Sample("7", 7.00, "oim_pf_7.jpg"),
    Sample("6.48", 6.48, "oim_pf_6.48.jpg"),
    Sample("6.02", 6.02, "oim_pf_6.02.jpg"),
    Sample("5.6", 5.60, "oim_pf_5.6.jpg"),
    Sample("5.25", 5.25, "oim_pf_5.25.jpg"),
    # oim_pf_5.jpg is the separately exported legend, not the pole figure.
    Sample("5", 5.00, "oim_pf_5_1.jpg"),
)


LEGEND_PATTERN = re.compile(
    r"^\s*([<>])\s*([0-9.]+)\s+Color rgb:\s*(\d+)\s+(\d+)\s+(\d+)\s*$",
    flags=re.MULTILINE,
)


def font(size: int, *, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = (
        Path("C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf"),
        Path("C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf"),
    )
    for candidate in candidates:
        if candidate.exists():
            return ImageFont.truetype(str(candidate), size=size)
    raise FileNotFoundError("Arial or Calibri was not found in C:/Windows/Fonts")


def parse_legend(path: Path) -> list[LegendBin]:
    text = path.read_text(encoding="utf-8", errors="replace")
    bins = [
        LegendBin(
            match.group(1),
            float(match.group(2)),
            (int(match.group(3)), int(match.group(4)), int(match.group(5))),
        )
        for match in LEGEND_PATTERN.finditer(text)
    ]
    if len(bins) != 8:
        raise ValueError(f"Expected 8 OIM color bins in {path}; found {len(bins)}")
    if bins[0].operator != "<" or bins[-1].operator != ">":
        raise ValueError(f"Unexpected lower/upper bin operators in {path}")
    if abs(bins[1].threshold - 1.0) > 1e-9:
        raise ValueError(f"Random-reference threshold 1.000 is missing in {path}")
    return bins


def centered_text(
    draw: ImageDraw.ImageDraw,
    center_x: int,
    y: int,
    text: str,
    text_font: ImageFont.FreeTypeFont,
    *,
    fill: str = "black",
) -> None:
    box = draw.textbbox((0, 0), text, font=text_font)
    draw.text((center_x - (box[2] - box[0]) / 2, y), text, font=text_font, fill=fill)


def centered_lines(
    draw: ImageDraw.ImageDraw,
    center_x: int,
    y: int,
    lines: list[str],
    text_font: ImageFont.FreeTypeFont,
    *,
    line_gap: int = 5,
    fill: str = "black",
) -> int:
    current_y = y
    for line in lines:
        centered_text(draw, center_x, current_y, line, text_font, fill=fill)
        box = draw.textbbox((0, 0), line, font=text_font)
        current_y += box[3] - box[1] + line_gap
    return current_y


def interval_label(bins: list[LegendBin], index: int) -> str:
    if index == 0:
        return f"<{bins[0].threshold:.3f}"
    if index == len(bins) - 1:
        return f">{bins[-1].threshold:.3f}"
    return f"{bins[index - 1].threshold:.3f}-{bins[index].threshold:.3f}"


def draw_intensity_legend(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    bins: list[LegendBin],
    legend_font: ImageFont.FreeTypeFont,
) -> None:
    centered_lines(
        draw,
        x + 78,
        y,
        ["OIM intensity", "1 = random", "auto bins"],
        legend_font,
        line_gap=2,
    )
    bar_x = x + 2
    bar_y = y + 100
    bar_width = 34
    block_height = 67
    for display_row, bin_index in enumerate(reversed(range(len(bins)))):
        entry = bins[bin_index]
        y0 = bar_y + display_row * block_height
        y1 = y0 + block_height
        draw.rectangle((bar_x, y0, bar_x + bar_width, y1), fill=entry.color, outline="black", width=1)
        label = interval_label(bins, bin_index)
        label_box = draw.textbbox((0, 0), label, font=legend_font)
        label_y = y0 + (block_height - (label_box[3] - label_box[1])) / 2 - label_box[1]
        draw.text((bar_x + bar_width + 8, label_y), label, font=legend_font, fill="black")


def prepare_panel(path: Path, target_size: int, *, expect_size: tuple[int, int]) -> Image.Image:
    image = Image.open(path).convert("RGB")
    if image.size != expect_size:
        image.close()
        raise ValueError(f"Unexpected size for {path}: expected {expect_size}, found {image.size}")
    if image.size == (target_size, target_size):
        return image
    resized = image.resize((target_size, target_size), Image.Resampling.LANCZOS)
    image.close()
    return resized


def build_matrix(source_root: Path, output_dir: Path) -> tuple[Path, Path]:
    ipf_dir = source_root / "oim_ipf"
    pf_dir = source_root / "oim_pf"
    odf_dir = source_root / "oim_odf"
    required = [ipf_dir / "tuli.jpg"]
    for sample in SAMPLES:
        required.extend(
            (
                ipf_dir / f"oim_ipf_{sample.stem}.jpg",
                pf_dir / sample.pf_filename,
                odf_dir / f"oim_odf_{sample.stem}.jpg",
                pf_dir / f"{sample.stem}.txt",
                odf_dir / f"{sample.stem}.txt",
            )
        )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required source files:\n" + "\n".join(missing))

    panel_size = 1024
    outer_margin = 48
    row_label_width = 285
    column_width = 1215
    column_gap = 36
    header_height = 390
    panel_label_height = 55
    row_gap = 28
    footer_height = 92
    row_height = panel_label_height + panel_size + row_gap
    canvas_width = (
        2 * outer_margin
        + row_label_width
        + 3 * column_width
        + 2 * column_gap
    )
    canvas_height = header_height + len(SAMPLES) * row_height + footer_height
    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")
    draw = ImageDraw.Draw(canvas)

    title_font = font(56, bold=True)
    column_font = font(44, bold=True)
    header_font = font(27)
    row_font = font(43, bold=True)
    row_subfont = font(31)
    panel_font = font(38, bold=True)
    legend_font = font(22)
    footer_font = font(27)

    data_x0 = outer_margin + row_label_width
    column_starts = [data_x0 + i * (column_width + column_gap) for i in range(3)]
    column_centers = [start + column_width // 2 for start in column_starts]

    centered_text(
        draw,
        canvas_width // 2,
        12,
        "OIM EBSD texture evolution: IPF, {0001} pole figure, and ODF",
        title_font,
    )
    headers = (
        ("IPF-[001] / A3", ["Color reference: sample A3 (map normal)", "Map horizontal: A2 = AD"]),
        ("{0001} PF (c-axis)", ["Top: +A1 | left: +A2 = AD", "Center: A3 (map normal)"]),
        ("ODF sections", ["Bunge: phi1 horizontal, PHI vertical", "phi2 = 0-60 deg; interval = 10 deg"]),
    )
    for index, (heading, sublines) in enumerate(headers):
        centered_text(draw, column_centers[index], 83, heading, column_font)
        centered_lines(draw, column_centers[index], 139, sublines, header_font, line_gap=3)

    # Add the common Ti-Hex IPF color key exported by OIM (triangle only).
    key_source = Image.open(ipf_dir / "tuli.jpg").convert("RGB")
    if key_source.size != (417, 433):
        key_source.close()
        raise ValueError(f"Unexpected IPF legend size: {key_source.size}")
    key_crop = key_source.crop((0, 228, 235, 386))
    key_source.close()
    key_display = key_crop.resize((390, 262), Image.Resampling.LANCZOS)
    key_crop.close()
    key_x = column_centers[0] - key_display.width // 2
    canvas.paste(key_display, (key_x, 192))
    key_display.close()

    # Compact PF coordinate guide.
    pf_center = (column_centers[1], 290)
    radius = 76
    draw.ellipse(
        (
            pf_center[0] - radius,
            pf_center[1] - radius,
            pf_center[0] + radius,
            pf_center[1] + radius,
        ),
        outline="black",
        width=3,
    )
    draw.line((pf_center[0], pf_center[1] - radius, pf_center[0], pf_center[1]), fill="black", width=2)
    draw.line((pf_center[0] - radius, pf_center[1], pf_center[0], pf_center[1]), fill="black", width=2)
    centered_text(draw, pf_center[0], pf_center[1] - radius - 34, "+A1", header_font)
    draw.text((pf_center[0] - radius - 92, pf_center[1] - 15), "+A2", font=header_font, fill="black")
    centered_text(draw, pf_center[0], pf_center[1] - 15, "A3", header_font)

    centered_lines(
        draw,
        column_centers[2],
        236,
        ["Seven sections per panel", "Harmonic: L=16, smoothing=5 deg", "Sample symmetry: Triclinic"],
        header_font,
        line_gap=8,
    )
    draw.line((outer_margin, header_height - 1, canvas_width - outer_margin, header_height - 1), fill="#808080", width=2)

    transformed_panels: list[tuple[Image.Image, int, int]] = []
    for row_index, sample in enumerate(SAMPLES):
        row_y = header_height + row_index * row_height
        image_y = row_y + panel_label_height
        row_center_y = image_y + panel_size // 2
        centered_text(draw, outer_margin + row_label_width // 2, row_center_y - 55, f"{sample.diameter_mm:.2f} mm", row_font)
        centered_text(
            draw,
            outer_margin + row_label_width // 2,
            row_center_y + 5,
            f"area reduction",
            row_subfont,
            fill="#333333",
        )
        centered_text(
            draw,
            outer_margin + row_label_width // 2,
            row_center_y + 45,
            f"{sample.area_reduction_percent:.2f}%",
            row_subfont,
            fill="#333333",
        )

        paths = (
            ipf_dir / f"oim_ipf_{sample.stem}.jpg",
            pf_dir / sample.pf_filename,
            odf_dir / f"oim_odf_{sample.stem}.jpg",
        )
        expected_sizes = ((1024, 1024), (1280, 1280), (1024, 1024))
        pf_bins = parse_legend(pf_dir / f"{sample.stem}.txt")
        odf_bins = parse_legend(odf_dir / f"{sample.stem}.txt")
        legends = (None, pf_bins, odf_bins)

        for column_index, (path, expected_size, legend_bins) in enumerate(
            zip(paths, expected_sizes, legends, strict=True)
        ):
            panel_index = row_index * 3 + column_index
            panel_x = column_starts[column_index]
            draw.text(
                (panel_x + 4, row_y + 4),
                f"({ascii_lowercase[panel_index]})",
                font=panel_font,
                fill="black",
            )
            panel = prepare_panel(path, panel_size, expect_size=expected_size)
            canvas.paste(panel, (panel_x, image_y))
            draw.rectangle(
                (panel_x - 1, image_y - 1, panel_x + panel_size, image_y + panel_size),
                outline="#666666",
                width=1,
            )
            transformed_panels.append((panel, panel_x, image_y))
            if legend_bins is not None:
                draw_intensity_legend(
                    draw,
                    panel_x + panel_size + 12,
                    image_y + 110,
                    legend_bins,
                    legend_font,
                )

        separator_y = row_y + row_height - row_gap // 2
        draw.line(
            (outer_margin, separator_y, canvas_width - outer_margin, separator_y),
            fill="#D0D0D0",
            width=1,
        )

    footer_y = header_height + len(SAMPLES) * row_height + 10
    centered_text(
        draw,
        canvas_width // 2,
        footer_y,
        "PF and ODF use specimen-specific OIM auto bins: compare peak positions directly, not identical colors across rows.",
        footer_font,
        fill="#333333",
    )

    # Verify every in-memory panel against its pasted pixels before lossless export.
    for panel, x, y in transformed_panels:
        pasted = canvas.crop((x, y, x + panel_size, y + panel_size))
        if ImageChops.difference(panel, pasted).getbbox() is not None:
            raise RuntimeError(f"Pixel verification failed at ({x}, {y})")
        panel.close()

    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / "oim_ipf_pf_odf_diameter_matrix.png"
    tiff_path = output_dir / "oim_ipf_pf_odf_diameter_matrix.tif"
    canvas.save(png_path, format="PNG", dpi=(600, 600), optimize=True)
    canvas.save(tiff_path, format="TIFF", dpi=(600, 600), compression="tiff_lzw")
    canvas.close()
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
        default=repository / "results" / "oim_ipf_pf_odf_matrix",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    png, tiff = build_matrix(args.source_root.resolve(), args.output_dir.resolve())
    print(f"Created: {png}")
    print(f"Created: {tiff}")
    print("Rows: 7.00, 6.48, 6.02, 5.60, 5.25, 5.00 mm")
    print("Columns: IPF-[001]/A3, {0001} PF, ODF phi2 sections")
    print("Verified: 18 panels paired correctly; no panel was rotated or cropped")
