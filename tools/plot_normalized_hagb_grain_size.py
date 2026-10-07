"""Plot a separate, geometrically normalized version of manuscript Fig. 4.

Uses existing per-grain MTEX output; does not rerun reconstruction or modify
source data, original plots, or manuscript figure links.
"""
from pathlib import Path
import csv
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'results/mtex_grain_size_distribution_matrix'
ORIGINAL = ROOT / 'results/mtex_fig6_component_gallery'
OUT = ROOT / 'results/hagb_grain_size_geometric_normalization_2026-09-27'
STATES = [('7d', 7.00), ('6.48d', 6.50), ('6.02d', 6.02),
          ('5.6d', 5.60), ('5.25d', 5.25), ('5d', 5.00)]
COLORS = [(0.40, 0.40, 0.40), (0.00, 0.45, 0.70),
          (0.34, 0.71, 0.91), (0.00, 0.62, 0.45),
          (0.90, 0.62, 0.00), (0.84, 0.37, 0.00)]


def read_csv(path):
    with path.open(newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    guarded = list(SOURCE.glob('hagb_grain_size*.csv'))
    guarded += list((ORIGINAL / 'montages').glob('hagb_grain_size_six_state_montage.*'))
    guarded += [ORIGINAL / '07_hagb_grain_size/hagb_grain_size_summary_used.csv']
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in guarded}
    grains = read_csv(SOURCE / 'hagb_grain_size_by_grain.csv')
    summary = {r['sample']: r for r in read_csv(
        ORIGINAL / '07_hagb_grain_size/hagb_grain_size_summary_used.csv')}
    raw_hist = read_csv(SOURCE / 'hagb_grain_size_histograms.csv')
    first_hist = sorted([r for r in raw_hist if r['sample'] == '7d'],
                        key=lambda r: float(r['bin_lower_um']))
    # Keep Fig. 4's common 2 um bins and x-range. Re-bin individual normalized
    # diameters, rather than simply changing labels or scaling mean values.
    edges = np.array([float(first_hist[0]['bin_lower_um'])] +
                     [float(r['bin_upper_um']) for r in first_hist])
    assert np.allclose(np.diff(edges), 2.0)
    centers = (edges[:-1] + edges[1:]) / 2
    xfit = np.linspace(0.05, edges[-1], 1000)
    states, derived_grains, derived_hist, derived_summary = [], [], [], []
    for sample, diameter in STATES:
        rows = [r for r in grains if r['sample'] == sample]
        original = summary[sample]
        observed = np.array([float(r['ecd_um']) for r in rows])
        assert len(rows) == int(original['grain_count'])
        assert all(r['variant'] == 'raw' and int(r['num_pixels']) >= 5 for r in rows)
        assert float(original['hagb_threshold_deg']) == 15
        assert np.isclose(observed.mean(), float(original['ecd_number_mean_um']), atol=1e-12)
        factor = np.sqrt(diameter / 7.0)
        normalized = observed * factor
        reduction = (1 - (diameter / 7.0) ** 2) * 100
        counts, _ = np.histogram(normalized, bins=edges)
        assert counts.sum() == len(rows)
        frequency = counts / len(rows) * 100
        assert np.isclose(frequency.sum(), 100)
        mu, sigma = np.log(normalized).mean(), np.log(normalized).std(ddof=0)
        assert np.isclose(mu, float(original['number_lognormal_mu']) + np.log(factor))
        assert np.isclose(sigma, float(original['number_lognormal_sigma']))
        fit = 100 * 2 * np.exp(-0.5 * ((np.log(xfit) - mu) / sigma) ** 2) / (xfit * sigma * np.sqrt(2 * np.pi))
        raw_counts, _ = np.histogram(observed, bins=edges)
        expected = sorted([r for r in raw_hist if r['sample'] == sample],
                          key=lambda r: float(r['bin_lower_um']))
        assert np.allclose(raw_counts / len(rows) * 100,
                           [float(r['number_frequency_percent']) for r in expected])
        metadata = dict(sample=sample, source_diameter_mm=float(original['diameter_mm']),
                        nominal_diameter_mm=diameter, nominal_reduction_percent=reduction,
                        geometric_scale_factor=factor)
        for row, value in zip(rows, normalized):
            derived_grains.append({**metadata, 'grain_id': int(row['grain_id']),
                'num_pixels': int(row['num_pixels']), 'boundary_touching': row['boundary_touching'],
                'observed_ecd_um': float(row['ecd_um']), 'normalized_ecd_um': value})
        for lower, upper, count, freq in zip(edges[:-1], edges[1:], counts, frequency):
            derived_hist.append({**metadata, 'bin_lower_um': lower, 'bin_upper_um': upper,
                                 'grain_count': int(count), 'number_frequency_percent': freq})
        derived_summary.append({**metadata, 'grain_count': len(rows),
            'hagb_threshold_deg': 15, 'min_grain_pixels': 5,
            'observed_number_mean_ecd_um': observed.mean(),
            'normalized_number_mean_ecd_um': normalized.mean(),
            'normalized_lognormal_mu': mu, 'lognormal_sigma': sigma})
        states.append(dict(diameter=diameter, reduction=reduction, mean=normalized.mean(),
                           frequency=frequency, fit=fit))
    # Match the manuscript table at its reported precision.
    assert [round(s['mean'], 2) for s in states] == [5.64, 6.02, 5.69, 5.51, 4.90, 4.83]
    ymax = max(5, 5 * np.ceil(1.12 * max(max(s['frequency'].max(), s['fit'].max()) for s in states) / 5))
    plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 9,
                         'axes.linewidth': 0.7, 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(2, 3, figsize=(11.2, 7.5))
    fig.subplots_adjust(left=.062, right=.987, bottom=.08, top=.86, wspace=.20, hspace=.29)
    fig.suptitle('Geometrically normalized HAGB grain ECD', y=.976, fontsize=13, fontweight='bold')
    fig.text(.5, .936, r'$d_{\mathrm{ECD,norm}} = d_{\mathrm{ECD,obs}}\sqrt{d/d_0}$;  '
             r'$d_0 = 7.00$ mm;  HAGB $\geq 15^\circ$;  minPixel = 5;  bin width = 2 $\mu$m',
             ha='center', fontsize=9)
    for index, (ax, state, color) in enumerate(zip(axes.flat, states, COLORS)):
        bars = ax.bar(centers, state['frequency'], width=1.68, color=color,
                      alpha=.74, edgecolor='.16', linewidth=.35, label='Histogram')
        line, = ax.plot(xfit, state['fit'], color='.08', linewidth=1.35,
                        label='Descriptive lognormal fit')
        if index == 0:
            ax.legend(handles=[bars, line], loc='upper right', fontsize=7, frameon=False)
        ax.text(.96, .72 if index == 0 else .91,
                f"Mean normalized ECD = {state['mean']:.2f} " + r'$\mu$m',
                transform=ax.transAxes, ha='right', va='top', fontsize=8,
                color=color, fontweight='bold')
        ax.set_title(f"({chr(97+index)}) {state['reduction']:.2f}% | {state['diameter']:.2f} mm",
                     fontsize=10, fontweight='bold', pad=7)
        ax.set(xlim=(0, edges[-1]), ylim=(0, ymax),
               xlabel=r'Normalized HAGB grain ECD ($\mu$m)', ylabel='Number frequency (%)')
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.tick_params(direction='out', top=True, right=True, labelsize=8, length=3)
    base = OUT / 'hagb_grain_size_geometrically_normalized_six_state'
    for ext in ('png', 'tif', 'svg'):
        options = {'pil_kwargs': {'compression': 'tiff_lzw'}} if ext == 'tif' else {}
        fig.savefig(base.with_suffix('.' + ext), dpi=600, facecolor='white', **options)
    plt.close(fig)
    write_csv(OUT / 'normalized_hagb_by_grain.csv', derived_grains)
    write_csv(OUT / 'normalized_hagb_histograms.csv', derived_hist)
    write_csv(OUT / 'normalized_hagb_summary.csv', derived_summary)
    assert all(digest(ROOT / path) == value for path, value in hashes.items())
    record = dict(source_sha256=hashes, originals_unchanged=True,
                  formula='ECD_norm = ECD_observed * sqrt(nominal_diameter_mm / 7.00)',
                  source_sample_to_nominal_diameter=dict(STATES),
                  bin_width_um=2, retained_grains=len(derived_grains),
                  assumptions=['uniform axisymmetric volume-preserving deformation',
                               'grain geometry follows macroscopic deformation'],
                  interpretation='Reference-geometry-normalized longitudinal-section ECD; not a measured 3D grain size.',
                  exclusions='No new exclusions; original >=5-pixel grain selection and boundary-touching grains retained.',
                  validation='All raw histograms and means reproduced; transformed lognormal parameters verified; normalized means match Table 3.')
    (OUT / 'provenance.json').write_text(json.dumps(record, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'figure': str(base.with_suffix('.png')), 'grain_count': len(derived_grains),
                      'normalized_means_um': [s['mean'] for s in states],
                      'source_files_unchanged': True}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
