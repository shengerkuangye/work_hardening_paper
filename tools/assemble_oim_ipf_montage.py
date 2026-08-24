"""Assemble the six OIM axial IPF maps into a publication-ready montage.

The source JPEGs are decoded and pasted at their native 1024 x 1024 pixel
size.  They are never rotated, cropped, resized, or overwritten.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw, ImageFont


@dataclass(frozen=True)
class Panel:
    filename: str
    diameter_mm: float
    letter: str

    @property
    def area_reduction_percent(self) -> float:
        initial_diameter_mm = 7.0
        return 100.0 * (1.0 - (self.diameter_mm / initial_diameter_mm) ** 2)


PANELS = (
    Panel("oim_ipf_7.jpg", 7.00, "a"),
    Panel("oim_ipf_6.48.jpg", 6.48, "b"),
    Panel("oim_ipf_6.02.jpg", 6.02, "c"),
    Panel("oim_ipf_5.6.jpg", 5.60, "d"),
    Panel("oim_ipf_5.25.jpg", 5.25, "e"),
    Panel("oim_ipf_5.jpg", 5.00, "f"),
)


def find_font(size: int, *, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = (
        Path("C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf"),
        Path("C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf"),
    )
    for candidate in candidates:
        if candidate.exists():
            return ImageFont.truetype(str(candidate), size=size)
    raise FileNotFoundError("Arial or Calibri font was not found in C:/Windows/Fonts")


def centered_text(
    draw: ImageDraw.ImageDraw,
    center_x: int,
    y: int,
    text: str,
    font: ImageFont.FreeTypeFont,
    fill: str = "black",
) -> None:
    box = draw.textbbox((0, 0), text, font=font)
    width = box[2] - box[0]
    draw.text((center_x - width / 2, y), text, font=font, fill=fill)


def build_montage(source_dir: Path, output_dir: Path) -> tuple[Path, Path]:
    missing = [panel.filename for panel in PANELS if not (source_dir / panel.filename).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing source maps: {', '.join(missing)}")

    loaded = [(panel, Image.open(source_dir / panel.filename).convert("RGB")) for panel in PANELS]
    sizes = {image.size for _, image in loaded}
    if len(sizes) != 1:
        raise ValueError(f"All source maps must have the same dimensions; found {sorted(sizes)}")

    image_width, image_height = next(iter(sizes))
    columns, rows = 3, 2
    outer_margin = 60
    column_gap = 32
    row_gap = 40
    top_header = 145
    panel_header = 68
    bottom_margin = 45

    canvas_width = 2 * outer_margin + columns * image_width + (columns - 1) * column_gap
    canvas_height = (
        top_header
        + rows * (panel_header + image_height)
        + (rows - 1) * row_gap
        + bottom_margin
    )
    montage = Image.new("RGB", (canvas_width, canvas_height), "white")
    draw = ImageDraw.Draw(montage)

    title_font = find_font(46, bold=True)
    direction_font = find_font(35)
    panel_font = find_font(42, bold=False)
    panel_letter_font = find_font(44, bold=True)

    centered_text(draw, canvas_width // 2, 13, "OIM IPF maps", title_font)

    arrow_y = 112
    arrow_start = canvas_width // 2 - 315
    arrow_end = canvas_width // 2 + 315
    centered_text(draw, canvas_width // 2, 65, "Axial direction (AD)", direction_font)
    draw.line((arrow_start, arrow_y, arrow_end, arrow_y), fill="black", width=5)
    draw.polygon(
        ((arrow_end, arrow_y), (arrow_end - 24, arrow_y - 13), (arrow_end - 24, arrow_y + 13)),
        fill="black",
    )

    paste_positions: list[tuple[Panel, Image.Image, int, int]] = []
    for index, (panel, image) in enumerate(loaded):
        row, column = divmod(index, columns)
        x = outer_margin + column * (image_width + column_gap)
        header_y = top_header + row * (panel_header + image_height + row_gap)
        image_y = header_y + panel_header

        letter = f"({panel.letter})"
        label = (
            f"{panel.diameter_mm:.2f} mm  |  "
            f"area reduction {panel.area_reduction_percent:.2f}%"
        )
        draw.text((x + 8, header_y + 9), letter, font=panel_letter_font, fill="black")
        letter_box = draw.textbbox((0, 0), letter, font=panel_letter_font)
        letter_width = letter_box[2] - letter_box[0]
        draw.text((x + 22 + letter_width, header_y + 11), label, font=panel_font, fill="black")

        montage.paste(image, (x, image_y))
        paste_positions.append((panel, image, x, image_y))

    # Assert pixel-for-pixel identity of every pasted map before export.
    for panel, source_image, x, y in paste_positions:
        pasted = montage.crop((x, y, x + image_width, y + image_height))
        if ImageChops.difference(source_image, pasted).getbbox() is not None:
            raise RuntimeError(f"Pixel verification failed for {panel.filename}")

    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / "oim_ipf_diameter_montage.png"
    tiff_path = output_dir / "oim_ipf_diameter_montage.tif"
    montage.save(png_path, format="PNG", dpi=(600, 600), optimize=True)
    montage.save(tiff_path, format="TIFF", dpi=(600, 600), compression="tiff_lzw")

    for _, image in loaded:
        image.close()
    return png_path, tiff_path


def parse_args() -> argparse.Namespace:
    repository = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-dir",
        type=Path,
        default=repository / "data" / "ebsd_kpl_250221_7_df" / "images_oim",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repository / "results" / "oim_ipf_montage",
    )
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    png, tiff = build_montage(arguments.source_dir.resolve(), arguments.output_dir.resolve())
    print(f"Created: {png}")
    print(f"Created: {tiff}")
    print("Panel order: 7.00, 6.48, 6.02, 5.60, 5.25, 5.00 mm")
    print("Verified: all six maps were pasted without rotation, resizing, or cropping")
