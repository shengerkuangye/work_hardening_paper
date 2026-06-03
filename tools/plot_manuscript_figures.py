from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
TENSILE = ROOT / "data" / "tensile_data" / "gr4b23271_cold_deformation_2_raw_csv"
TENSILE_ROOT = ROOT / "data" / "tensile_data"


COLORS = {
    "blue": "#1f77b4",
    "red": "#d62728",
    "green": "#2ca02c",
    "orange": "#ff7f0e",
    "purple": "#9467bd",
    "brown": "#8c564b",
    "gray": "#444444",
    "grid": "#d9d9d9",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def fmt(value: float, digits: int = 2) -> str:
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


class Svg:
    def __init__(self, width: int, height: int) -> None:
        self.width = width
        self.height = height
        self.parts: list[str] = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            "<style>",
            "text{font-family:Arial,Helvetica,sans-serif;fill:#222;font-size:16px}",
            ".small{font-size:13px}.axis{stroke:#222;stroke-width:1.4}.grid{stroke:#d9d9d9;stroke-width:1}",
            ".line{fill:none;stroke-width:2.6;stroke-linejoin:round;stroke-linecap:round}",
            ".marker{stroke:#fff;stroke-width:1.3}",
            "</style>",
        ]

    def line(self, x1: float, y1: float, x2: float, y2: float, color: str, cls: str = "", width: float | None = None) -> None:
        extra = f' class="{cls}"' if cls else ""
        sw = f' stroke-width="{width}"' if width else ""
        self.parts.append(f'<line x1="{fmt(x1)}" y1="{fmt(y1)}" x2="{fmt(x2)}" y2="{fmt(y2)}" stroke="{color}"{sw}{extra}/>')

    def text(self, x: float, y: float, text: str, anchor: str = "start", cls: str = "", rotate: int | None = None) -> None:
        extra = f' class="{cls}"' if cls else ""
        transform = f' transform="rotate({rotate} {fmt(x)} {fmt(y)})"' if rotate is not None else ""
        self.parts.append(f'<text x="{fmt(x)}" y="{fmt(y)}" text-anchor="{anchor}"{extra}{transform}>{text}</text>')

    def polyline(self, points: list[tuple[float, float]], color: str, cls: str = "line") -> None:
        pts = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in points)
        self.parts.append(f'<polyline points="{pts}" stroke="{color}" class="{cls}"/>')

    def circle(self, x: float, y: float, r: float, color: str) -> None:
        self.parts.append(f'<circle cx="{fmt(x)}" cy="{fmt(y)}" r="{fmt(r)}" fill="{color}" class="marker"/>')

    def rect(self, x: float, y: float, w: float, h: float, fill: str = "none", stroke: str = "#222") -> None:
        self.parts.append(f'<rect x="{fmt(x)}" y="{fmt(y)}" width="{fmt(w)}" height="{fmt(h)}" fill="{fill}" stroke="{stroke}" stroke-width="1"/>')

    def save(self, path: Path) -> None:
        self.parts.append("</svg>")
        path.write_text("\n".join(self.parts), encoding="utf-8")


def plot_axes(svg: Svg, x0: int, y0: int, w: int, h: int, xlim: tuple[float, float], ylim: tuple[float, float], xticks: list[float], yticks: list[float], xlabel: str, ylabel: str) -> tuple:
    xmin, xmax = xlim
    ymin, ymax = ylim

    def sx(x: float) -> float:
        return x0 + (x - xmin) / (xmax - xmin) * w

    def sy(y: float) -> float:
        return y0 + h - (y - ymin) / (ymax - ymin) * h

    for t in xticks:
        x = sx(t)
        svg.line(x, y0, x, y0 + h, COLORS["grid"], "grid")
        svg.text(x, y0 + h + 22, fmt(t, 0), "middle", "small")
    for t in yticks:
        y = sy(t)
        svg.line(x0, y, x0 + w, y, COLORS["grid"], "grid")
        svg.text(x0 - 10, y + 5, fmt(t, 0), "end", "small")
    svg.line(x0, y0 + h, x0 + w, y0 + h, "#222", "axis")
    svg.line(x0, y0, x0, y0 + h, "#222", "axis")
    svg.text(x0 + w / 2, y0 + h + 50, xlabel, "middle")
    svg.text(x0 - 55, y0 + h / 2, ylabel, "middle", rotate=-90)
    return sx, sy


def add_legend(svg: Svg, entries: list[tuple[str, str]], x: int, y: int) -> None:
    for i, (label, color) in enumerate(entries):
        yy = y + i * 24
        svg.line(x, yy - 4, x + 28, yy - 4, color, width=2.8)
        svg.circle(x + 14, yy - 4, 4, color)
        svg.text(x + 36, yy, label, "start", "small")


def make_mechanical_summary() -> None:
    rows = read_rows(TENSILE_ROOT / "gr4b23271_lab_tensile_by_diameter.csv")
    x = [float(r["cold_reduction_percent_reference_vs_7mm"]) for r in rows]
    rp = [float(r["Rp0.2_MPa_mean"]) for r in rows]
    rm = [float(r["Rm_MPa_mean"]) for r in rows]
    elong = [float(r["total_elongation_A_percent_mean"]) for r in rows]
    z = [float(r["reduction_of_area_Z_percent_mean"]) for r in rows]

    svg = Svg(980, 720)
    svg.text(90, 38, "(a) Strength", "start")
    sx, sy = plot_axes(svg, 90, 65, 810, 245, (0, 50), (500, 1100), [0, 10, 20, 30, 40, 50], [500, 700, 900, 1100], "Cold reduction (%)", "Strength (MPa)")
    for values, color in [(rp, COLORS["blue"]), (rm, COLORS["red"])]:
        pts = [(sx(a), sy(b)) for a, b in zip(x, values)]
        svg.polyline(pts, color)
        for px, py in pts:
            svg.circle(px, py, 5, color)
    add_legend(svg, [("Rp0.2", COLORS["blue"]), ("Rm", COLORS["red"])], 735, 88)

    svg.text(90, 392, "(b) Ductility", "start")
    sx2, sy2 = plot_axes(svg, 90, 420, 810, 220, (0, 50), (0, 50), [0, 10, 20, 30, 40, 50], [0, 10, 20, 30, 40, 50], "Cold reduction (%)", "Ductility (%)")
    for values, color in [(elong, COLORS["green"]), (z, COLORS["orange"])]:
        pts = [(sx2(a), sy2(b)) for a, b in zip(x, values)]
        svg.polyline(pts, color)
        for px, py in pts:
            svg.circle(px, py, 5, color)
    add_legend(svg, [("A", COLORS["green"]), ("Z", COLORS["orange"])], 735, 443)

    svg.save(FIGURES / "gr4b23271_mechanical_property_summary.svg")


def moving_slope(x: list[float], y: list[float], window: int) -> list[tuple[float, float]]:
    half = window // 2
    out: list[tuple[float, float]] = []
    for i in range(half, len(x) - half):
        xs = x[i - half : i + half + 1]
        ys = y[i - half : i + half + 1]
        xm = sum(xs) / len(xs)
        ym = sum(ys) / len(ys)
        denom = sum((v - xm) ** 2 for v in xs)
        if denom <= 0:
            continue
        slope = sum((a - xm) * (b - ym) for a, b in zip(xs, ys)) / denom
        out.append((x[i], slope))
    return out


def resample_linear(x: list[float], y: list[float], step: float) -> tuple[list[float], list[float]]:
    if len(x) < 2:
        return [], []
    x0 = math.ceil(x[0] / step) * step
    x1 = math.floor(x[-1] / step) * step
    grid: list[float] = []
    values: list[float] = []
    j = 0
    n = int(round((x1 - x0) / step)) + 1
    for i in range(max(0, n)):
        gx = x0 + i * step
        while j + 1 < len(x) and x[j + 1] < gx:
            j += 1
        if j + 1 >= len(x):
            break
        x_left, x_right = x[j], x[j + 1]
        y_left, y_right = y[j], y[j + 1]
        if x_right == x_left:
            continue
        ratio = (gx - x_left) / (x_right - x_left)
        grid.append(gx)
        values.append(y_left + ratio * (y_right - y_left))
    return grid, values


def moving_average(values: list[float], window: int) -> list[float]:
    if not values:
        return []
    window = max(3, window | 1)
    half = window // 2
    out: list[float] = []
    running = 0.0
    counts = 0
    queue: list[float] = []
    for i, value in enumerate(values):
        queue.append(value)
        running += value
        counts += 1
        if len(queue) > window:
            running -= queue.pop(0)
            counts -= 1
        if i < half:
            out.append(sum(values[: i + half + 1]) / len(values[: i + half + 1]))
        elif i >= len(values) - half:
            tail = values[i - half :]
            out.append(sum(tail) / len(tail))
        else:
            out.append(running / counts)
    return out


def central_difference(x: list[float], y: list[float], half_span: int) -> list[tuple[float, float]]:
    out: list[tuple[float, float]] = []
    if len(x) <= 2 * half_span:
        return out
    for i in range(half_span, len(x) - half_span):
        dx = x[i + half_span] - x[i - half_span]
        if dx <= 0:
            continue
        out.append((x[i], (y[i + half_span] - y[i - half_span]) / dx))
    return out


def solve_linear_system(a: list[list[float]], b: list[float]) -> list[float]:
    n = len(b)
    matrix = [row[:] + [rhs] for row, rhs in zip(a, b)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(matrix[r][col]))
        if abs(matrix[pivot][col]) < 1e-12:
            raise ValueError("Singular matrix")
        matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
        pivot_value = matrix[col][col]
        for c in range(col, n + 1):
            matrix[col][c] /= pivot_value
        for r in range(n):
            if r == col:
                continue
            factor = matrix[r][col]
            for c in range(col, n + 1):
                matrix[r][c] -= factor * matrix[col][c]
    return [matrix[i][n] for i in range(n)]


def polynomial_fit(x: list[float], y: list[float], degree: int) -> list[float]:
    degree = min(degree, len(x) - 1)
    size = degree + 1
    powers = [sum(value ** p for value in x) for p in range(2 * degree + 1)]
    lhs = [[powers[i + j] for j in range(size)] for i in range(size)]
    rhs = [sum((xv ** i) * yv for xv, yv in zip(x, y)) for i in range(size)]
    return solve_linear_system(lhs, rhs)


def polynomial_value(coeffs: list[float], x: float) -> float:
    return sum(c * (x ** i) for i, c in enumerate(coeffs))


def polynomial_derivative(coeffs: list[float], x: float) -> float:
    return sum(i * c * (x ** (i - 1)) for i, c in enumerate(coeffs) if i > 0)


def display_sample_name(sample: str, diameter: str) -> str:
    if "M" in sample:
        return "TA4-M"
    value = float(diameter)
    if abs(value - 6.48) < 0.05:
        return "TA4-Y-6.5"
    return f"TA4-Y-{fmt(value, 2)}"


def make_hardening_rate() -> None:
    reps = read_rows(FIGURES / "gr4b23271_representative_tensile_curves.csv")
    summary_rows = {r["sample"]: r for r in read_rows(TENSILE / "gr4b23271_cold_deformation_tensile_summary.csv")}
    series: list[tuple[str, float, list[tuple[float, float]]]] = []
    derived_rows: list[dict[str, str]] = []

    for rep in reps:
        sample = rep["sample"]
        raw = read_rows(TENSILE / rep["source_file"])
        meta = summary_rows[sample]
        rp_eng = float(meta["strain_at_Rp0.2"]) if meta["strain_at_Rp0.2"] else 0.0
        rp_true = math.log1p(max(rp_eng, 0.0))
        fracture_true_strain = float(meta["fracture_true_strain_last_percent"]) / 100.0
        xs: list[float] = []
        ys: list[float] = []
        last_x = -1.0
        for r in raw:
            try:
                x = float(r["Ture Strain"])
                y = float(r["Ture Stress"])
            except (KeyError, ValueError):
                continue
            if x < rp_true + 0.001 or x > fracture_true_strain - 0.001 or x <= last_x:
                continue
            xs.append(x)
            ys.append(y)
            last_x = x
        if len(xs) < 80:
            continue
        grid_x, grid_y = resample_linear(xs, ys, 0.00025)
        if len(grid_x) < 12:
            continue
        smooth_window = max(7, min(41, len(grid_x) // 5 | 1, int(round(0.0025 / 0.00025)) | 1))
        smoothed_y = moving_average(grid_y, smooth_window)
        half_span = max(4, min(14, len(grid_x) // 8, int(round(0.002 / 0.00025))))
        slopes = central_difference(grid_x, smoothed_y, half_span)
        step = max(1, len(slopes) // 450)
        reduced = [(xval, max(theta, 0.0)) for xval, theta in slopes[::step]]
        label = display_sample_name(sample, rep["diameter_mm"])
        cold = float(rep["cold_reduction_percent_reference"])
        series.append((label, cold, reduced))
        for xval, theta in reduced:
            derived_rows.append(
                {
                    "sample": sample,
                    "cold_reduction_percent": f"{cold:.6g}",
                    "true_strain": f"{xval:.8g}",
                    "work_hardening_rate_MPa": f"{theta:.8g}",
                }
            )

    with (FIGURES / "gr4b23271_work_hardening_rate_representative_curves.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "cold_reduction_percent", "true_strain", "work_hardening_rate_MPa"])
        writer.writeheader()
        writer.writerows(derived_rows)

    xmax = max(max(p[0] for p in pts) for _, _, pts in series)
    ymin, ymax = 0.0, 12000.0
    svg = Svg(980, 620)
    x_axis_max = min(0.38, math.ceil(xmax * 20) / 20)
    xticks = [0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35]
    sx, sy = plot_axes(svg, 95, 70, 780, 455, (0, x_axis_max), (ymin, ymax), xticks, [0, 3000, 6000, 9000, 12000], "True strain", "Work hardening rate (MPa)")
    palette = [COLORS["gray"], COLORS["blue"], COLORS["red"], COLORS["green"], COLORS["orange"], COLORS["purple"]]
    legend_entries: list[tuple[str, str]] = []
    for (label, cold, pts), color in zip(series, palette):
        plot_pts = [(sx(x), sy(theta)) for x, theta in pts if 0 <= x <= x_axis_max and ymin <= theta <= ymax]
        if len(plot_pts) > 1:
            svg.polyline(plot_pts, color)
        legend_entries.append((f"{label} ({fmt(cold, 1)}%)", color))
    add_legend(svg, legend_entries, 665, 96)
    svg.save(FIGURES / "gr4b23271_work_hardening_rate_representative_curves.svg")


def main() -> None:
    FIGURES.mkdir(exist_ok=True)
    make_mechanical_summary()
    make_hardening_rate()


if __name__ == "__main__":
    main()
