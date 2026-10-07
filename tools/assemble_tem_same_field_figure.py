#!/usr/bin/env python3
"""Assemble existing experimental frames; never synthesize or alter contrast.

Foil-edge markers establish field overlap, not dislocation identity.
Diffraction labels denote provisional row families, not verified operating g.
"""
from pathlib import Path
import argparse
import base64
import hashlib
import io
import json
import subprocess
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'data/9_DF/DF'
OUT = ROOT / 'results/tem_same_field_2026-10-01'
AUDIT = ROOT / 'results/tem_field_correspondence_2026-10-01'
RECORDS = json.loads((AUDIT / 'feature_correspondence_audit.json').read_text())
TRIPLES = [('0026', '0028', '0019'), ('0029', '0031', '0020'), ('0032', '0034', '0022')]


def matrix(a, b):
    r = next(r for r in RECORDS if r['source'] == a and r['target'] == b)
    m = np.eye(3)
    m[:2] = r['matrix']
    return m


def source(aid, dp=False):
    return next(SRC.glob(aid + ('*Camera Ceta.tif' if dp else '*STEM BF-S.tif')))


def embedded(path):
    with Image.open(path) as im:
        buffer = io.BytesIO()
        im.convert('RGB').save(buffer, format='PNG')
    return 'data:image/png;base64,' + base64.b64encode(buffer.getvalue()).decode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state', choices=['5', '7'], default='5')
    args = parser.parse_args()
    is_m7 = args.state == '7'
    out = ROOT / 'results/tem_m7_multicondition_2026-10-01' if is_m7 else OUT
    triples = [('0009', '0011', '0003'), ('0012', '0014', '0004'), ('0015', '0017', '0005')] if is_m7 else TRIPLES
    out.mkdir(parents=True, exist_ok=True)
    side, gap, margin = 1000, 24, 24
    width, height = 3096, 2900
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>']
    manifest = []
    for col, (low, high, dp) in enumerate(triples):
        x = margin + col * (side + gap)
        H = matrix(low, high)
        G = None if is_m7 else (np.eye(3) if low == '0026' else matrix('0026', low))
        for row, (aid, is_dp, y, size) in enumerate([(low, False, 62, side), (high, False, 1138, side), (dp, True, 2208, 660)]):
            xx = x if not is_dp else x + 170
            path = source(aid, is_dp)
            label = chr(97 + row * 3 + col)
            desc = f'({label}) {aid}'
            if is_dp:
                desc = f'({label}) DP {aid}'
            svg.append(f'<text x="{x}" y="{y-15}" font-family="DejaVu Sans" font-size="45" fill="black">{desc}</text>')
            svg.append(f'<image x="{xx}" y="{y}" width="{size}" height="{size}" href="{embedded(path)}"/>')
            if is_dp:
                family = ['(0002)', '{10<tspan text-decoration="overline">1</tspan>0}', '{10<tspan text-decoration="overline">1</tspan>1}'][col]
                svg.append(f'<text x="{x+side}" y="{y-15}" text-anchor="end" font-family="DejaVu Sans" font-size="45">{family}*</text>')
            if row == 0:
                corners = np.array([[0,0,1], [2048,0,1], [2048,2048,1], [0,2048,1]]) @ np.linalg.inv(H).T
                points = ' '.join(f'{x+p[0]*side/2048:.2f},{y+p[1]*side/2048:.2f}' for p in corners)
                svg.append(f'<polygon points="{points}" stroke="#ffdc44" stroke-width="5" fill="none"/>')
                # Only the 5 mm series has sufficiently supported cross-condition
                # landmark mapping. Never transfer its marker to the 7 mm series.
                if G is not None:
                    point = G @ np.array([1343.5, 1918.6, 1])
                    px, py = x+point[0]*side/2048, y+point[1]*side/2048
                    svg.append(f'<circle cx="{px}" cy="{py}" r="17" fill="none" stroke="#00dcff" stroke-width="5"/>')
                    svg.append(f'<text x="{px+24}" y="{py-12}" font-family="DejaVu Sans" font-size="44" fill="#00dcff" stroke="black" stroke-width="2" paint-order="stroke">L</text>')
            manifest.append(dict(panel=label, acquisition=aid, source=str(path.relative_to(ROOT)), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    svg.append('</svg>')
    path = out / ('fig_m7_multicondition_diffraction.svg' if is_m7 else 'fig12_y5_same_field_diffraction.svg')
    path.write_text('\n'.join(svg))
    subprocess.run(['/usr/bin/inkscape', str(path), '--export-type=png', '--export-width=3600', '--export-filename='+str(path.with_suffix('.png'))], check=True)
    marker_note = 'Yellow polygon: measured coverage of high-magnification image in the same column.'
    if not is_m7:
        marker_note = 'L: shared foil-edge landmark; ' + marker_note
    field_status = 'Within-condition low/high maps verified. Cross-condition registration and individual dislocation identity unconfirmed; no common-feature labels.' if is_m7 else 'Foil-edge landmarks support overlapping fields; individual dislocation identity unconfirmed.'
    (out / 'provenance.json').write_text(json.dumps(dict(state='M-7.00' if is_m7 else 'Y-5.00', panels=manifest, image_processing='Full original frames, proportional resizing only; no warping, contrast edits or synthetic defects.', markers=marker_note, field_status=field_status, diffraction='Provisional systematic-row family labels (*); no assertion of actual selected g or single-dislocation extinction.', source_annotation_mapping=triples), ensure_ascii=False, indent=2)+'\n')
    print(path.with_suffix('.png'))


if __name__ == '__main__':
    main()
