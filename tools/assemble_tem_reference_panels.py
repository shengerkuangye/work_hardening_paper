#!/usr/bin/env python3
"""Assemble traceable TEM resource pairs; never modify experimental inputs.

Display operations: RGB conversion and proportional LANCZOS reduction only.
No crop, intensity adjustment, denoising, rotation, registration or synthesis.
Source diffraction associations are annotations, NOT validated g conditions.
"""
import csv
import hashlib
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'data/9_DF/DF'
OUT = ROOT / 'results/tem_9df_reference_panels'
FONT = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
STATES = {'7MM': 'M-7.00', '5MM': 'Y-5.00'}
GROUPS = [('7MM', 'C1', 1, (6, 7, 8)), ('7MM', 'C2', 3, (9, 10, 11)),
          ('7MM', 'C3', 4, (12, 13, 14)), ('7MM', 'C4', 5, (15, 16, 17)),
          ('5MM', 'C1', 18, (23, 24, 25)), ('5MM', 'C2', 19, (26, 27, 28)),
          ('5MM', 'C3', 20, (29, 30, 31)), ('5MM', 'C4', 22, (32, 33, 34))]
SPECS = [('A', 'a', 9, 'BF-S', 'overview'),
         ('A', 'b', 11, 'BF-S', 'curved/intersecting line contrast'),
         ('A', 'c', 11, 'DF-S', 'same-acquisition complementary detector'),
         ('A', 'd', 26, 'BF-S', 'overview'),
         ('A', 'e', 28, 'BF-S', 'curved/intersecting line contrast'),
         ('A', 'f', 28, 'DF-S', 'same-acquisition complementary detector'),
         ('B', 'a', 11, 'BF-S', 'condition C2; diffraction tilt mismatch'),
         ('B', 'b', 17, 'BF-S', 'condition C4; g unindexed'),
         ('B', 'c', 28, 'BF-S', 'condition C2; g unindexed'),
         ('B', 'd', 34, 'BF-S', 'condition C4; g unindexed')]


def find(aid, suffix):
    paths = list(SRC.glob(f'{aid:04d} - *{suffix}'))
    assert len(paths) == 1, (aid, suffix, paths)
    return paths[0]


def meta(path):
    with Image.open(path) as im:
        root = ET.fromstring(im.tag_v2[34683])
        dims = im.size
    return dict(path=str(path.relative_to(ROOT)), width=dims[0], height=dims[1],
                alpha_deg=float(root.findtext('StageSettings/StagePosition/Tilt/Alpha'))*180/math.pi,
                beta_deg=float(root.findtext('StageSettings/StagePosition/Tilt/Beta'))*180/math.pi,
                time_utc=root.findtext('Acquisition/AcquisitionStartDatetime'),
                detector_xml=root.findtext('BinaryResult/Detector'),
                pixel_size=float(root.findtext('BinaryResult/PixelSize/X')))


def csv_out(name, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with (OUT / name).open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def label(canvas, xy, text, size=24, fill='#20252b'):
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype(FONT, size)
    assert xy[0] + draw.textlength(text, font=font) <= canvas.width - 8, text
    draw.text(xy, text, font=font, fill=fill)


def paste(canvas, path, xy, side):
    with Image.open(path) as source:
        im = source.convert('RGB')
    im.thumbnail((side, side), Image.Resampling.LANCZOS)
    canvas.paste(im, xy)


def export(canvas, name, tiff=False):
    canvas.save(OUT / (name + '.png'), dpi=(300, 300))
    if tiff:
        canvas.save(OUT / (name + '.tif'), compression='tiff_lzw', dpi=(300, 300))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'selected_full_frames').mkdir(exist_ok=True)
    paired, channels, lookup = [], [], {}
    for group, cond, did, aids in GROUPS:
        dm = meta(find(did, 'Camera Ceta.tif'))
        for aid in aids:
            sm = meta(find(aid, 'STEM BF-S.tif'))
            da, db = sm['alpha_deg'] - dm['alpha_deg'], sm['beta_deg'] - dm['beta_deg']
            # 0.05 degrees is a bookkeeping flag, not a physical two-beam tolerance.
            flag = 'tilt_record_close_only' if max(abs(da), abs(db)) < .05 else 'tilt_record_mismatch'
            row = dict(state=STATES[group], source_group=group, section='longitudinal',
                       condition_label=cond, stem_id=f'{aid:04d}', diffraction_id=f'{did:04d}',
                       bf_tiff=sm['path'], diffraction_tiff=dm['path'],
                       annotated_diffraction_png=str(find(did, 'Ceta-0.png').relative_to(ROOT)),
                       pair_basis='original diffraction PNG acquisition-number annotation',
                       stem_alpha_deg=sm['alpha_deg'], stem_beta_deg=sm['beta_deg'],
                       diffraction_alpha_deg=dm['alpha_deg'], diffraction_beta_deg=dm['beta_deg'],
                       delta_alpha_deg=da, delta_beta_deg=db, tilt_record_flag=flag,
                       g_status='unindexed; actual STEM excitation condition unverified',
                       field_identity='gross landmarks compatible with overlap; individual dislocations not tracked',
                       stem_time_utc=sm['time_utc'], diffraction_time_utc=dm['time_utc'])
            paired.append(row)
            lookup[aid] = row
            for ch in ['BF-S', 'DF-S', 'DF-O', 'HAADF']:
                p = find(aid, f'STEM {ch}.tif')
                cm = meta(p)
                channels.append(dict(state=STATES[group], stem_id=f'{aid:04d}', condition_label=cond,
                                     diffraction_id=f'{did:04d}', filename_channel=ch,
                                     tiff=str(p.relative_to(ROOT)), mrc=str(p.with_suffix('.mrc').relative_to(ROOT)),
                                     emd=str(find(aid, 'STEM.emd').relative_to(ROOT)),
                                     detector_xml=cm['detector_xml'],
                                     detector_note='filename/XML identity differs; preserve source labels' if ch == 'DF-O' else '',
                                     alpha_deg=cm['alpha_deg'], beta_deg=cm['beta_deg'],
                                     time_utc=cm['time_utc']))
    csv_out('diffraction_stem_pair_audit.csv', paired)
    csv_out('channel_file_pairs.csv', channels)

    selected = []
    source_hashes = {}
    for fig, panel, aid, channel, purpose in SPECS:
        path = find(aid, f'STEM {channel}.tif')
        row = lookup[aid]
        selected.append(dict(figure=fig, panel=panel, state=row['state'], stem_id=f'{aid:04d}',
                             channel=channel, condition_label=row['condition_label'],
                             source_tiff=str(path.relative_to(ROOT)), diffraction_id=row['diffraction_id'],
                             purpose=purpose, processing='full exported frame; proportional resizing only',
                             interpretation_limit='morphology only' if fig == 'A' else 'candidate comparison; no assigned g or b'))
        source_hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        # Full-resolution lossless viewing copy for formats not displayed by the UI.
        png = OUT / 'selected_full_frames' / f"{row['state']}_{aid:04d}_{channel}.png"
        with Image.open(path) as src:
            src.convert('RGB').save(png)
    csv_out('selected_panels.csv', selected)

    # A: morphology, two states x overview / high-magnification BF / paired DF.
    side, gap = 1000, 20
    canvas = Image.new('RGB', (3080, 2270), 'white')
    label(canvas, (20, 10), 'A | Morphology selection | TA4 | Longitudinal sections', 34)
    label(canvas, (20, 57), 'Overview BF-S  /  detail BF-S  /  same-acquisition DF-S. No density ranking.', 26)
    for i, (_, panel, aid, ch, _) in enumerate(SPECS[:6]):
        x, y = 20 + (i % 3)*(side+gap), 106 + (i // 3)*1080
        row = lookup[aid]
        label(canvas, (x, y), f"({panel}) {row['state']} | {aid:04d} | {ch}", 27)
        paste(canvas, find(aid, f'STEM {ch}.tif'), (x, y+42), side)
    label(canvas, (20, 2240), 'Full exported frames; original scale bars retained. Selection draft, not a quantitative comparison.', 23)
    export(canvas, 'A_morphology_selection', tiff=True)

    # B: same-channel, same-nominal-magnification comparison. Diffraction placed
    # beside, not on top of, each full real-space frame; association caveat visible.
    canvas = Image.new('RGB', (2440, 2090), 'white')
    label(canvas, (20, 12), 'B | Candidate dislocation-identification pairs | BF-S, 45k', 33)
    label(canvas, (20, 58), 'g UNINDEXED. Diffraction panels are source-linked records, not verified STEM conditions.', 24)
    for i, (_, panel, aid, ch, _) in enumerate(SPECS[6:]):
        x, y = 20 + (i % 2)*1210, 110 + (i // 2)*960
        row = lookup[aid]
        label(canvas, (x, y), f"({panel}) {row['state']} | {aid:04d} | {row['condition_label']}", 28)
        paste(canvas, find(aid, f'STEM {ch}.tif'), (x, y+48), 880)
        dx, dy = x+896, y+48
        paste(canvas, find(int(row['diffraction_id']), 'Ceta-0.png'), (dx, dy), 280)
        label(canvas, (dx, dy+292), f"DP {row['diffraction_id']}", 23)
        label(canvas, (dx, dy+326), 'Source annotation', 20)
        label(canvas, (dx, dy+368), 'Tilt: STEM / DP', 20)
        label(canvas, (dx, dy+401), f"a: {row['stem_alpha_deg']:.2f} / {row['diffraction_alpha_deg']:.2f}", 19)
        label(canvas, (dx, dy+431), f"b: {row['stem_beta_deg']:.2f} / {row['diffraction_beta_deg']:.2f}", 19)
        flag = row['tilt_record_flag']
        label(canvas, (dx, dy+479), 'TILT MISMATCH' if flag.endswith('mismatch') else 'Tilt records close', 20,
              '#ad4717' if flag.endswith('mismatch') else '#20252b')
        label(canvas, (dx, dy+518), 'g: not assigned', 21)
    label(canvas, (20, 2055), 'No individual-dislocation extinction claim; no slip-plane assignment. a/b beside DP denote stage tilt angles (degrees).', 22)
    export(canvas, 'B_diffraction_pair_candidates', tiff=True)

    # Complete review pairing: one row per source annotation, four columns.
    for group in STATES:
        canvas = Image.new('RGB', (2090, 2460), 'white')
        label(canvas, (18, 10), f"{STATES[group]} | Full resource pairing | Longitudinal | g unindexed", 27)
        label(canvas, (18, 49), 'Source diffraction PNG / BF-S 22.5k / BF-S 32k / BF-S 45k', 24)
        for r, (_, cond, did, aids) in enumerate(g for g in GROUPS if g[0] == group):
            y = 93 + r*587
            da = max(abs(lookup[a]['delta_alpha_deg']) for a in aids)
            db = max(abs(lookup[a]['delta_beta_deg']) for a in aids)
            label(canvas, (18, y), f"{cond} | DP {did:04d} -> {aids[0]:04d}-{aids[-1]:04d} | max |delta tilt|: a={da:.3f}, b={db:.3f} deg", 22)
            paths = [find(did, 'Ceta-0.png')] + [find(a, 'STEM BF-S.tif') for a in aids]
            for c, path in enumerate(paths):
                paste(canvas, path, (18+c*518, y+34), 500)
                label(canvas, (18+c*518, y+538), path.name.split(' - ')[0] + (' | Diffraction' if c == 0 else ' | BF-S'), 22)
        export(canvas, f'{STATES[group]}_complete_pairing')

    for rel, expected in source_hashes.items():
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == expected
    assert len(paired) == 24 and len(channels) == 96 and len(selected) == 10
    for row in channels:
        assert all((ROOT / row[k]).is_file() for k in ['tiff', 'mrc', 'emd'])
    for path in OUT.glob('*.png'):
        with Image.open(path) as im:
            im.verify()
    report = dict(paired_stem_acquisitions=len(paired), channel_records=len(channels),
                  selected_panel_records=len(selected),
                  tilt_close_records=sum(r['tilt_record_flag']=='tilt_record_close_only' for r in paired),
                  tilt_mismatch_records=sum(r['tilt_record_flag']=='tilt_record_mismatch' for r in paired),
                  tilt_audit_threshold_deg=.05,
                  threshold_note='bookkeeping only; closeness is not proof of identical excitation or field',
                  selected_source_sha256_unchanged=source_hashes,
                  manuscript_modified=False, raw_source_modified=False)
    (OUT / 'build_validation.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k != 'selected_source_sha256_unchanged'}, indent=2))


if __name__ == '__main__':
    main()
