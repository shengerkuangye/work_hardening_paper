#!/usr/bin/env python3
"""Vector annotations over unmodified experimental TIFF exports.

This is an analysis/selection figure, not a confirmed Burgers-vector map.
All feature coordinates use the original 2048 x 2048 STEM frame.
No registration, contrast enhancement, tracing, denoising or synthesis.
"""
from pathlib import Path
import base64
import csv
import hashlib
import io
import json
import math
import html
from PIL import Image, ImageDraw, ImageFont

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'data/9_DF/DF'
OUT=ROOT/'results/tem_extinction_annotated_2026-09-30'
MEASURE=ROOT/'results/tem_extinction_selection_2026-09-30/systematic_row_measurements.json'
BLUE='#00c5ef'
YELLOW='#ffe062'
PANELS=[
    dict(panel='a',state='M-7.00',aid='0014',dp='0004',family='10-10',features=[
        ('A1',(820,905),(605,830))]),
    dict(panel='b',state='M-7.00',aid='0011',dp='0003',family='0002',features=[
        ('B1',(1270,644),(1110,510)),('B2',(778,1120),(570,1030))]),
    dict(panel='c',state='M-7.00',aid='0017',dp='0005',family='10-11',features=[
        ('C1',(1460,1050),(1240,940))]),
    dict(panel='d',state='Y-5.00',aid='0031',dp='0020',family='10-10',features=[
        ('D1',(360,1790),(150,1650))]),
    dict(panel='e',state='Y-5.00',aid='0028',dp='0019',family='0002',features=[
        ('E1',(410,536),(225,440)),('E2',(621,652),(400,735))]),
    dict(panel='f',state='Y-5.00',aid='0034',dp='0022',family='10-11',features=[
        ('F1',(1590,1390),(1370,1280))]),
]


def source(aid,kind):
    return next(SRC.glob(aid+('*STEM BF-S.tif' if kind=='stem' else '*Camera Ceta.tif')))


def png_uri(path):
    with Image.open(path) as im:
        b=io.BytesIO();im.convert('RGB').save(b,format='PNG')
    return 'data:image/png;base64,'+base64.b64encode(b.getvalue()).decode()


def txt(x,y,s,size=24,fill='#17202a',weight='normal',anchor='start',halo=False):
    return (f'<text x="{x}" y="{y}" font-family="DejaVu Sans, sans-serif" '
            f'font-size="{size}" font-weight="{weight}" fill="{fill}" text-anchor="{anchor}" '
            + ('stroke="#18212a" stroke-width="5" paint-order="stroke" stroke-linejoin="round" ' if halo else '')
            + '>'+html.escape(s)+'</text>')


def family_text(x,y,family,size=34):
    # An actual overbar, not a minus sign standing in for a crystallographic bar.
    if family=='0002':return txt(x,y,'(0002)',size)
    return (f'<text x="{x}" y="{y}" font-family="DejaVu Sans" font-size="{size}" fill="#17202a">'
            '{10<tspan text-decoration="overline">1</tspan>'+family[-1]+'}</text>')


def arrow(x1,y1,x2,y2,color,width=4):
    vx,vy=x2-x1,y2-y1;n=math.hypot(vx,vy);ux,uy=vx/n,vy/n
    # Filled triangular head ends precisely at the selected source coordinate.
    bx,by=x2-17*ux,y2-17*uy
    return (f'<path d="M{x1},{y1} L{bx},{by}" stroke="#17202a" stroke-width="{width+3}" fill="none"/>'
            f'<path d="M{x1},{y1} L{bx},{by}" stroke="{color}" stroke-width="{width}" fill="none"/>'
            f'<polygon points="{x2},{y2} {bx-7*uy},{by+7*ux} {bx+7*uy},{by-7*ux}" fill="{color}" stroke="#17202a" stroke-width="1"/>')


def inspect_targets():
    font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',22)
    sheet=Image.new('RGB',(1200,6*460),'white');draw=ImageDraw.Draw(sheet)
    for row,p in enumerate(PANELS):
        with Image.open(source(p['aid'],'stem')) as im:
            for col,(label,(x,y),_) in enumerate(p['features']):
                box=(x-200,y-200,x+200,y+200)
                crop=im.convert('RGB').crop(box)
                # Only a perimeter box at the target, retaining underlying pixels.
                cd=ImageDraw.Draw(crop);cd.rectangle((186,186,214,214),outline='red',width=2)
                sheet.paste(crop,(col*600+10,row*460+44))
                draw.text((col*600+10,row*460+8),f"{label} | {p['aid']} | source ({x},{y})",font=font,fill='black')
    sheet.save(OUT/'annotation_target_audit.png')


def build(annotated=True):
    rows={r['diffraction_id']:r for r in json.loads(MEASURE.read_text())}
    margin,gap,side=36,24,1000
    width=2*margin+3*side+2*gap
    rowheight=1400
    height=2*rowheight+120
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
           '<rect width="100%" height="100%" fill="white"/>',
           txt(margin,48,'TA4 | Dislocation-contrast comparison',34,weight='bold'),
           txt(margin,84,'Working annotation  /  longitudinal section  /  BF-STEM  /  full original frames',25)]
    for i,p in enumerate(PANELS):
        x=margin+(i%3)*(side+gap);y=120+(i//3)*rowheight
        parts += [txt(x,y+32,f"({p['panel']})  {p['state']}  |  {p['aid']}",30,weight='bold'),
                  family_text(x+650,y+32,p['family'],30),txt(x+805,y+32,'row',27)]
        iy=y+55
        parts.append(f'<image x="{x}" y="{iy}" width="{side}" height="{side}" xlink:href="{png_uri(source(p["aid"],"stem"))}"/>')
        parts.append(f'<rect x="{x}" y="{iy}" width="{side}" height="{side}" fill="none" stroke="#b9bec4" stroke-width="1"/>')
        color=BLUE if p['family']=='0002' else YELLOW
        if annotated:
            for label,(tx,ty),(sx,sy) in p['features']:
                sc=side/2048
                parts += [arrow(x+sx*sc,iy+sy*sc,x+tx*sc,iy+ty*sc,color),
                          txt(x+sx*sc-8,iy+sy*sc-14,label,31,fill=color,weight='bold',halo=True)]
        dy=iy+side+18;ds=290
        parts.append(f'<image x="{x}" y="{dy}" width="{ds}" height="{ds}" xlink:href="{png_uri(source(p["dp"],"dp"))}"/>')
        # Show the measured row direction alongside, not an assumed direct beam
        # or a verified STEM excitation g. No hkil indices are painted on spots.
        fit=rows[p['dp']];vx,vy=fit['step_x'],fit['step_y'];n=math.hypot(vx,vy)
        parts.append(arrow(x+ds+54,dy+208,x+ds+54+85*vx/n,dy+208+85*vy/n,'#2774a6',3))
        parts += [txt(x+ds+30,dy+36,'DP '+p['dp']+' | provisional row indexing',24,weight='bold'),
                  txt(x+ds+30,dy+74,f"Measured d = {fit['d_nm']:.3f} nm",24),
                  txt(x+ds+30,dy+112,'Arrow: systematic-row direction',22)]
    parts.append('</svg>')
    name='fig12_tem_annotated_review' if annotated else 'fig12_tem_clean_layout'
    (OUT/(name+'.svg')).write_text('\n'.join(parts),encoding='utf-8')


def manifest():
    rows=[]
    for p in PANELS:
        for label,target,start in p['features']:
            rows.append(dict(panel=p['panel'],image_id=p['aid'],diffraction_id=p['dp'],
                             feature_label=label,source_x=target[0],source_y=target[1],
                             label_x=start[0],label_y=start[1],
                             role='retained-line/c-component candidate under provisional indexing' if p['family']=='0002' else 'comparison line feature',
                             same_dislocation_cross_panel='not established',
                             burgers_vector='unassigned'))
    with (OUT/'annotation_coordinates.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    paths=[source(p['aid'],'stem') for p in PANELS]+[source(p['dp'],'dp') for p in PANELS]
    record=dict(purpose='User-requested annotation and assembly of existing experimental images',
                c_plus_a_status='candidate for further analysis; no unique Burgers-vector assignment',
                source_pixels='full frames, embedded lossless PNG conversions; no local enhancement or synthetic reconstruction',
                annotation_semantics='Panel-specific labels; not a completed extinction/visibility matrix',
                source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (OUT/'provenance.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')


if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=True)
    inspect_targets();build(True);build(False);manifest()
    print(OUT)
