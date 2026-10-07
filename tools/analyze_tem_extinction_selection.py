#!/usr/bin/env python3
"""Read raw diffraction MRC arrays and measure systematic-row spacings.

Run with PYTHONPATH=.codex_tmp/tem_analysis_deps python3 ... .
Source data are read only; labels are analysis labels, not final indexing.
"""
from pathlib import Path
import csv
import json
import struct
import hashlib
import xml.etree.ElementTree as ET
import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'data/9_DF/DF'
OUT = ROOT / 'results/tem_extinction_selection_2026-09-30'
FONT = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
IDS = ['0001', '0003', '0004', '0005', '0018', '0019', '0020', '0022']
# Seed points inspected in the raw diffraction peak audit; these select rows,
# not Burgers vectors. The central seed is only an order origin for spacing.
ROWS = {
    '0003': ((2041,2155),(100,357),'0002'),
    '0004': ((1861,2169),(323,-121),'10-10'),
    '0005': ((1997,2127),(386,44),'10-11'),
    '0019': ((1981,2058),(-68,364),'0002'),
    '0020': ((1955,2087),(340,57),'10-10'),
    '0022': ((1966,2103),(314,224),'10-11'),
}


def raw_mrc(path):
    with path.open('rb') as f:
        header = f.read(1024)
        nx, ny, nz, mode = struct.unpack_from('<4i', header)
        extra = struct.unpack_from('<i', header, 92)[0]
        f.seek(1024 + extra)
        data = np.fromfile(f, dtype={1: '<i2', 2: '<f4', 6: '<u2'}[mode], count=nx*ny*nz)
    return data.reshape(nz, ny, nx)[0].astype(float)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    all_peaks = {}
    for aid in IDS:
        raw = raw_mrc(next(SRC.glob(aid+'*Ceta.mrc')))
        smooth = gaussian_filter(raw, 2)
        ys, xs = np.where((smooth == maximum_filter(smooth, 41)) & (smooth > max(20, smooth.max()*.002)))
        keep = np.argsort(smooth[ys, xs])[::-1][:80]
        xs, ys = xs[keep], ys[keep]
        with Image.open(next(SRC.glob(aid+'*Ceta.tif'))) as im:
            root = ET.fromstring(im.tag_v2[34683])
            cal = float(root.findtext('BinaryResult/PixelSize/X')) / 1e9
            frame = im.convert('RGB').resize((1200, 1200), Image.Resampling.LANCZOS)
        points = []
        for x,y in zip(xs,ys):
            # Local background-subtracted centroid; detection on original MRC.
            patch = raw[y-10:y+11,x-10:x+11]
            if patch.shape != (21,21): continue
            weights = np.maximum(patch - np.median(patch),0)
            yy,xx=np.indices(patch.shape)
            cx=x-10+(weights*xx).sum()/weights.sum()
            cy=y-10+(weights*yy).sum()/weights.sum()
            points.append([float(cx),float(cy),float(smooth[y,x])])
        draw=ImageDraw.Draw(frame)
        font=ImageFont.truetype(FONT,17)
        for i,(x,y,v) in enumerate(points[:20]):
            x,y=x*1200/raw.shape[1],y*1200/raw.shape[0]
            draw.ellipse((x-9,y-9,x+9,y+9),outline='red',width=1)
            draw.text((x+10,y-12),str(i),font=font,fill='yellow')
        frame.save(OUT / (aid+'_peak_audit.png'))
        all_peaks[aid]={'calibration_nm_inv_per_pixel':cal, 'points_xy_intensity':points}
        print(aid, 'cal',cal,'top',np.round(points[:9],2).tolist())
    (OUT/'diffraction_peaks.json').write_text(json.dumps(all_peaks,indent=2)+'\n')


def fit_rows():
    all_peaks=json.loads((OUT/'diffraction_peaks.json').read_text())
    summaries=[]
    theory={'0002':2/.468,'10-10':np.sqrt(4/3)/.295,
            '10-11':np.sqrt(4/(3*.295**2)+1/.468**2)}
    for aid,(origin,step,family) in ROWS.items():
        pts=np.array(all_peaks[aid]['points_xy_intensity'])
        origin,step=np.array(origin),np.array(step)
        n=np.round((pts[:,:2]-origin)@step/(step@step)).astype(int)
        residual=np.linalg.norm(pts[:,:2]-origin-n[:,None]*step,axis=1)
        # Keep highest-intensity detection within each expected lattice order.
        ix=[]
        for order in sorted(set(n[residual<35])):
            options=np.where((n==order)&(residual<35))[0]
            ix.append(options[np.argmax(pts[options,2])])
        n=n[ix];p=pts[ix,:2]
        fit=np.linalg.lstsq(np.column_stack([np.ones(len(n)),n]),p,rcond=None)[0]
        rms=np.sqrt(np.mean(np.sum((p-n[:,None]*fit[1]-fit[0])**2,axis=1)))
        measured=np.linalg.norm(fit[1])*all_peaks[aid]['calibration_nm_inv_per_pixel']
        row={'diffraction_id':aid,'candidate_family':family,'n_spots':len(n),
             'step_pixels':float(np.linalg.norm(fit[1])), 'g_nm_inv':float(measured),
             'd_nm':float(1/measured),'theory_g_nm_inv':float(theory[family]),
             'relative_difference_percent':float((measured/theory[family]-1)*100),
             'fit_rms_pixels':float(rms),'origin_x':float(fit[0,0]),'origin_y':float(fit[0,1]),
             'step_x':float(fit[1,0]),'step_y':float(fit[1,1])}
        summaries.append(row)
        print(row)
    with (OUT/'systematic_row_measurements.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summaries[0]));w.writeheader();w.writerows(summaries)
    (OUT/'systematic_row_measurements.json').write_text(json.dumps(summaries,indent=2)+'\n')


def zone_check():
    peaks=json.loads((OUT/'diffraction_peaks.json').read_text())
    results=[]
    for aid,origin,basis in [('0001',(1972,2099),((325,-105),(56,176))),
                              ('0018',(1988,2093),((338,59),(-33,183)))]:
        p=np.array(peaks[aid]['points_xy_intensity'])[:,:2]
        basis=np.array(basis,dtype=float)
        hk=np.round((p-origin)@np.linalg.inv(basis))
        residual=np.linalg.norm(p-origin-hk@basis,axis=1)
        mask=residual<18
        fit=np.linalg.lstsq(np.column_stack([np.ones(mask.sum()),hk[mask]]),p[mask],rcond=None)[0]
        a,c=fit[1],2*fit[2]
        cal=peaks[aid]['calibration_nm_inv_per_pixel']
        results.append({'diffraction_id':aid,'fit_spots':int(mask.sum()),
                        'a_star_row_g_nm_inv':float(np.linalg.norm(a)*cal),
                        '0002_g_nm_inv':float(np.linalg.norm(c)*cal),
                        'angle_deg':float(np.degrees(np.arccos(a@c/np.linalg.norm(a)/np.linalg.norm(c)))),
                        'row_ratio_c_to_a':float(np.linalg.norm(c)/np.linalg.norm(a)),
                        'basis_a_pixels':a.tolist(),'basis_c_pixels':c.tolist(),
                        'interpretation':'Near <1-210>-type zone metric; representative basal reciprocal axis, not unique signed indices.'})
    (OUT/'zone_metric_check.json').write_text(json.dumps(results,indent=2)+'\n')
    print('ZONE',results)


def selection_panels():
    font=ImageFont.truetype(FONT,24)
    small=ImageFont.truetype(FONT,20)
    specs=[('M-7.00',['0014','0011','0017'],['0004','0003','0005']),
           ('Y-5.00',['0031','0028','0034'],['0020','0019','0022'])]
    family=['{10-10} row','(0002) row','{10-11} row']
    # All panels use full, unaltered TIFF exports; only proportional resizing.
    for state,ids,dps in specs:
        w,h=2190,1150
        canvas=Image.new('RGB',(w,h),'white');d=ImageDraw.Draw(canvas)
        d.text((18,12),state+' | extinction-based selection | 45k BF-S | longitudinal',font=font,fill='black')
        d.text((18,48),'Provisional row indexing from raw MRC spacings + zone-pattern metric; not a Burgers-vector assignment.',font=small,fill='black')
        for j,(aid,dp) in enumerate(zip(ids,dps)):
            x=18+725*j
            d.text((x,86),aid+' | '+family[j]+' | DP '+dp,font=font,fill='black')
            with Image.open(next(SRC.glob(aid+'*STEM BF-S.tif'))) as im:
                im.save(OUT/(aid+'_BF-S_full.png'))
                frame=im.convert('RGB').resize((700,700),Image.Resampling.LANCZOS)
            canvas.paste(frame,(x,125))
            with Image.open(next(SRC.glob(dp+'*Ceta.tif'))) as im:
                frame=im.convert('RGB').resize((285,285),Image.Resampling.LANCZOS)
            canvas.paste(frame,(x,845))
            d.text((x+302,855),'Diffraction source',font=small,fill='black')
            d.text((x+302,885),'Full frames; resized only',font=small,fill='black')
        canvas.save(OUT/(state+'_extinction_selection.png'))
    # Lower magnification full fields for manual correspondence, no registration
    # or local image enhancement has been applied.
    for state,ids in [('M-7.00',['0012','0009','0015']),('Y-5.00',['0029','0026','0032'])]:
        canvas=Image.new('RGB',(2190,840),'white');d=ImageDraw.Draw(canvas)
        d.text((18,12),state+' | matching low-magnification fields | 22.5k BF-S',font=font,fill='black')
        for j,aid in enumerate(ids):
            x=18+725*j
            d.text((x,54),aid+' | '+family[j],font=font,fill='black')
            with Image.open(next(SRC.glob(aid+'*STEM BF-S.tif'))) as im:
                frame=im.convert('RGB').resize((700,700),Image.Resampling.LANCZOS)
            canvas.paste(frame,(x,100))
        canvas.save(OUT/(state+'_overlapping_fields.png'))


def provenance():
    paths=[]
    for aid in IDS:
        paths.extend(SRC.glob(aid+'*Ceta.mrc'))
        paths.extend(SRC.glob(aid+'*Ceta.tif'))
    for aid in ['0009','0011','0012','0014','0015','0017','0026','0028','0029','0031','0032','0034']:
        paths.extend(SRC.glob(aid+'*STEM BF-S.tif'))
    record={
        'purpose':'Select existing TEM/STEM images using provisional diffraction-row indexing; no new measurements requested.',
        'method':'MRC local maxima, manually selected systematic rows, integer-order least-squares spacing fit; independent two-dimensional zone-pattern metric check.',
        'nominal_lattice_nm':{'a':.295,'c':.468},
        'lattice_reference':'https://www.sciencedirect.com/science/article/pii/S0167577X11015072',
        'primary_pairs':{'M-7.00':['0014 / DP0004','0011 / DP0003'],
                         'Y-5.00':['0031 / DP0020','0028 / DP0019']},
        'third_conditions':{'M-7.00':'0017 / DP0005','Y-5.00':'0034 / DP0022'},
        'interpretation':[
            '0003 and 0019 match the 0002 row; 0004 and 0020 match {10-10}; 0005 and 0022 match {10-11}.',
            'These are provisional row-family identifications from recorded diffraction patterns, not certifications of actual STEM excitation error or strict two-beam conditions.',
            'Equivalent indices/signs are not uniquely assigned. Curly braces indicate reflection family, not a confirmed active slip plane.',
            '0002 is a useful primary g.b test for a-type vs c-component contrast; additional constraints are needed to uniquely distinguish c from c+a.',
            'Visual low-magnification landmarks support overlapping fields. Automated high-magnification SIFT matches were insufficient and rejected; no image registration, warping or same-dislocation labels were applied.',
            'TIFF stage-angle differences remain as recorded in the earlier pairing audit; these do not change the measured diffraction-row spacings.',
            'No completed per-dislocation visibility matrix or a-to-c+a population transition is claimed.'
        ],
        'image_processing':'Full TIFF frames; proportional resizing and external text only. Peak-audit images have explicit analysis overlays. No synthetic structures or selective enhancement.',
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(paths))}
    }
    (OUT/'selection_provenance.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')


if __name__ == '__main__':
    import sys
    if '--fit-only' not in sys.argv:
        main()
    fit_rows()
    zone_check()
    selection_panels()
    provenance()
