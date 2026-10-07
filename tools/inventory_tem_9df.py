#!/usr/bin/env python3
"""Read-only source inventory and unenhanced contact sheets for TEM batch 9_DF.

TIFF metadata come from FEI XML tag 34683. Contact sheets preserve the full
exported image, including original scale bars; only proportional resizing is
applied. File count and acquisition ID count are NOT independent specimen n.
"""
import csv
import json
import math
import re
import struct
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/9_DF/DF'
OUT = ROOT / 'results/tem_9df_inventory'
FONT = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
STATES = {'7MM': 'M-7.00', '5MM': 'Y-5.00'}  # Author confirmed, 2026-09-26.


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def metadata(path):
    with Image.open(path) as im:
        row = {'width_px': im.width, 'height_px': im.height,
               'export_mode': im.mode}
        xml = im.tag_v2.get(34683)
        if not xml:
            row['metadata_status'] = 'no_FEI_XML'
            return row
        root = ET.fromstring(xml)
        fields = {
            'instrument': 'Instrument/InstrumentModel',
            'voltage_V': 'Optics/AccelerationVoltage',
            'acquisition_start_utc': 'Acquisition/AcquisitionStartDatetime',
            'operating_mode': 'Optics/OperatingMode',
            'projector_mode': 'Optics/ProjectorMode',
            'camera_length_m': 'Optics/CameraLength',
            'detector': 'BinaryResult/Detector',
            'pixel_size_x': 'BinaryResult/PixelSize/X',
            'pixel_size_y': 'BinaryResult/PixelSize/Y',
        }
        row.update({key: root.findtext(xp, '') for key, xp in fields.items()})
        pixel = root.find('BinaryResult/PixelSize/X')
        row['pixel_unit'] = pixel.get('unit', '') if pixel is not None else ''
        for axis in ('X', 'Y'):
            row[f'field_width_{axis}_um'] = (
                float(row[f'pixel_size_{axis.lower()}'])
                * row['width_px' if axis == 'X' else 'height_px'] * 1e6
                if row['pixel_unit'] == 'm' else '')
        row['metadata_status'] = 'ok'
        return row


def sheet(paths, name, title, cols=3, side=500):
    padding, label = 18, 44
    width = cols * (side + padding) + padding
    height = 74 + math.ceil(len(paths) / cols) * (side + label + padding)
    canvas = Image.new('RGB', (width, height), 'white')
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype(FONT, 19)
    draw.text((padding, 12), title, font=font, fill='black')
    draw.text((padding, 38), 'Full exported frames; resized only; original scale bars retained.', font=font, fill='black')
    for i, path in enumerate(paths):
        x = padding + (i % cols) * (side + padding)
        y = 74 + (i // cols) * (side + label + padding)
        with Image.open(path) as source:
            frame = source.convert('RGB')
        frame.thumbnail((side, side), Image.Resampling.LANCZOS)
        canvas.paste(frame, (x, y))
        short = path.stem.replace(' - ', ' | ').replace('Camera Ceta', 'Ceta')
        short = short.replace('22500 x STEM ', '22.5k | ').replace('32000 x STEM ', '32k | ').replace('45000 x STEM ', '45k | ')
        for source_group, state in STATES.items():
            short = short.replace(source_group, state)
        draw.text((x, y + side + 5), short, font=font, fill='black')
    canvas.save(OUT / name)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    file_rows, tif_rows, acqs = [], [], defaultdict(list)
    for path in sorted(SOURCE.iterdir()):
        if not path.is_file():
            continue
        match = re.match(r'(\d{4}) - (\d+MM) (.+)', path.stem)
        row = {'path': str(path.relative_to(ROOT)), 'size_bytes': path.stat().st_size,
               'extension': path.suffix, 'acquisition_id': match[1] if match else '',
               'source_group': match[2] if match else '',
               'filename_description': match[3] if match else path.stem}
        file_rows.append(row)
        row['manuscript_state'] = STATES.get(row['source_group'], '')
        row['section'] = 'longitudinal'
        row['sample_mapping_basis'] = 'author confirmation, 2026-09-26'
        acqs[(row['source_group'], row['acquisition_id'])].append(path)
        if path.suffix == '.tif':
            tif_rows.append({**row, **metadata(path)})
        elif path.suffix == '.mrc':
            with path.open('rb') as handle:
                header = handle.read(1024)
            nx, ny, nz, mode = struct.unpack_from('<4i', header)
            row.update(mrc_nx=nx, mrc_ny=ny, mrc_nz=nz, mrc_mode=mode)
    write_csv(OUT / 'file_manifest.csv', file_rows)
    write_csv(OUT / 'tiff_metadata.csv', tif_rows)
    acq_rows = []
    for (group, aid), paths in sorted(acqs.items()):
        tifs = [r for r in tif_rows if r['source_group'] == group and r['acquisition_id'] == aid]
        acq_rows.append({
            'source_group': group, 'acquisition_id': aid,
            'manuscript_state': STATES[group], 'section': 'longitudinal',
            'acquisition_kind': 'STEM' if any('STEM' in p.name for p in paths) else 'Camera',
            'file_count': len(paths),
            'detectors': ';'.join(sorted({r['detector'] for r in tifs})),
            'acquisition_start_utc': ';'.join(sorted({r['acquisition_start_utc'] for r in tifs})),
            'tiff_paths': ';'.join(r['path'] for r in tifs),
            'interpretation': 'acquisition identifier, not independent specimen or field count',
        })
    write_csv(OUT / 'acquisition_manifest.csv', acq_rows)
    summary = {
        'file_count': len(file_rows),
        'bytes': sum(r['size_bytes'] for r in file_rows),
        'extensions': dict(Counter(r['extension'] for r in file_rows)),
        'source_groups': dict(Counter(r['source_group'] for r in file_rows)),
        'acquisitions_by_group_and_kind': dict(Counter(r['source_group'] + '_' + r['acquisition_kind'] for r in acq_rows)),
        'instrument_models': sorted({r['instrument'] for r in tif_rows}),
        'acceleration_voltages_V': sorted({r['voltage_V'] for r in tif_rows}),
        'acquisition_start_dates_utc': sorted({r['acquisition_start_utc'][:10] for r in tif_rows}),
        'detector_counts': dict(Counter(r['detector'] for r in tif_rows)),
        'pixel_sizes_by_group_detector_unit': sorted({(r['source_group'],r['detector'],r['pixel_size_x'],r['pixel_unit']) for r in tif_rows}),
        'missing_integer_ids_in_range': sorted(set(range(1,35)) - {int(r['acquisition_id']) for r in file_rows}),
        'confirmed_state_mapping': STATES,
        'confirmed_section': 'longitudinal',
        'confirmation_basis': 'author confirmation, 2026-09-26',
        'interpretation': 'Missing ID does not establish a missing experiment. Source-to-state mapping and longitudinal section confirmed; independent specimen/field repeats remain unknown.',
    }
    (OUT / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    for group in ('7MM', '5MM'):
        for channel in ('BF-S', 'DF-S', 'DF-O', 'HAADF'):
            paths = sorted(SOURCE.glob(f'*{group}*STEM {channel}.tif'))
            sheet(paths, f'{group}_{channel}_contact.png', f'{STATES[group]} ({group}) | longitudinal | STEM {channel} | {len(paths)} acquisitions')
    sheet(sorted(SOURCE.glob('*Ceta-0.png')), 'diffraction_annotated_contact.png', 'Source-annotated diffraction exports', cols=4)
    sheet(sorted(SOURCE.glob('*Camera Ceta.tif')), 'camera_contact.png', 'Camera TIFF exports (includes diffraction and one real-space image)', cols=3)
    # Source PNG labels transcribed after visual inspection, not crystallographic indices.
    links = [('7MM','0001','0006;0007;0008','6-8'),
             ('7MM','0003','0009;0010;0011','9-11'),
             ('7MM','0004','0012;0013;0014','12-14'),
             ('7MM','0005','0015;0016;0017','15-17'),
             ('5MM','0018','0023;0024;0025','23-25'),
             ('5MM','0019','0026;0027;0028','26-28'),
             ('5MM','0020','0029;0030;0031','29-31'),
             ('5MM','0022','0032;0033;0034','32-34')]
    write_csv(OUT / 'diffraction_image_links.csv', [
        dict(source_group=g, manuscript_state=STATES[g], diffraction_acquisition_id=d, stem_acquisition_ids=a,
             original_png_label=l, status='source annotation; field identity and diffraction indexing not independently verified')
        for g,d,a,l in links])
    panels = ['0011 - 7MM 45000 x STEM BF-S.tif',
              '0011 - 7MM 45000 x STEM DF-S.tif',
              '0003 - 7MM 520 mm Camera Ceta-0.png',
              '0028 - 5MM 45000 x STEM BF-S.tif',
              '0028 - 5MM 45000 x STEM DF-S.tif',
              '0019 - 5MM 520 mm Camera Ceta-0.png']
    sheet([SOURCE / p for p in panels], 'tem_source_group_comparison_draft.png',
          'Review draft | M-7.00 (top) / Y-5.00 (bottom) | longitudinal | BF-S / DF-S / diffraction', cols=3, side=650)
    write_csv(OUT / 'comparison_panels.csv', [
        dict(panel=chr(97+i), path=str((SOURCE/p).relative_to(ROOT)),
             selection_reason='visible line contrast; illustrative only, not a density ranking',
             processing='full exported frame; proportional resizing only')
        for i,p in enumerate(panels)])
    print(json.dumps({k:v for k,v in summary.items() if k != 'pixel_sizes_by_group_detector_unit'}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
