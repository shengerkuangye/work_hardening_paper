#!/usr/bin/env python3
"""Common-field review using measured magnification maps and edge landmarks.

All outputs are original-intensity crops; no images are geometrically warped.
Cross-condition maps are edge-anchored estimates, not per-dislocation matches.
"""
from pathlib import Path
import json
import numpy as np
from PIL import Image,ImageDraw,ImageFont

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'data/9_DF/DF'
OUT=ROOT/'results/tem_field_correspondence_2026-10-01'
FONT='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
records=json.loads((OUT/'feature_correspondence_audit.json').read_text())


def matrix(a,b):
    r=next(r for r in records if r['source']==a and r['target']==b)
    m=np.eye(3);m[:2]=r['matrix'];return m


def original(aid):
    return Image.open(next(SRC.glob(aid+'*STEM BF-S.tif'))).convert('RGB')


def main():
    font=ImageFont.truetype(FONT,23);small=ImageFont.truetype(FONT,18)
    H26=matrix('0026','0028');H29=matrix('0029','0031');H32=matrix('0032','0034')
    G29=matrix('0026','0029');G32=matrix('0026','0032')
    high=[np.eye(3),H29@G29@np.linalg.inv(H26),H32@G32@np.linalg.inv(H26)]
    ids=['0028','0031','0034'];frames=[original(aid) for aid in ids]
    output=[]
    for row,y in enumerate([800,1200,1650]):
        canvas=Image.new('RGB',(1280,1510),'white');d=ImageDraw.Draw(canvas)
        d.text((16,12),'Predicted common regions | original pixels, no registration warp',font=font,fill='black')
        d.text((16,45),'Positions use foil-edge geometry; these are not confirmed same-dislocation matches.',font=small,fill='black')
        for k,x in enumerate([750,1100,1500]):
            for j,(aid,M,im) in enumerate(zip(ids,high,frames)):
                cx,cy=(M@np.array([x,y,1]))[:2]
                cx,cy=round(cx),round(cy)
                box=(cx-200,cy-200,cx+200,cy+200)
                xx=16+j*422;yy=100+k*455
                d.text((xx,yy),f'{aid} | R{row*3+k+1} | ({cx},{cy})',font=small,fill='black')
                crop=im.crop(box)
                canvas.paste(crop,(xx,yy+32))
                output.append(dict(image_id=aid,region=f'R{row*3+k+1}',box=list(box),center=[cx,cy],
                                   nominal_reference_image='0028',nominal_reference_point=[x,y],
                                   position_basis='edge-landmark similarity fit composed with measured within-condition magnification maps'))
        canvas.save(OUT/f'common_roi_review_{row+1}.png')
    (OUT/'common_roi_coordinates.json').write_text(json.dumps(output,indent=2)+'\n')
    (OUT/'high_magnification_estimated_maps.json').write_text(json.dumps({aid:M.tolist() for aid,M in zip(ids,high)},indent=2)+'\n')
    # Full low-magnification frames: high-mag footprint and the same three
    # stable foil-edge landmarks. These landmarks are not dislocations.
    lowids=['0026','0029','0032'];Gs=[np.eye(3),G29,G32];Hs=[H26,H29,H32]
    canvas=Image.new('RGB',(2490,1030),'white');d=ImageDraw.Draw(canvas)
    d.text((20,16),'Y-5.00 | same-field localization before dislocation analysis',font=font,fill='black')
    landmarks=[('L1',(1198.6,1743.4)),('L2',(1343.5,1918.6)),('L3',(1554.1,1882.4))]
    for j,(aid,G,H) in enumerate(zip(lowids,Gs,Hs)):
        x=20+j*825;y=95;scale=800/2048
        d.text((x,59),f'{aid} | 22.5k BF-S | high-mag {ids[j]}',font=font,fill='black')
        canvas.paste(original(aid).resize((800,800),Image.Resampling.LANCZOS),(x,y))
        corners=np.array([[0,0,1],[2048,0,1],[2048,2048,1],[0,2048,1]])@np.linalg.inv(H).T
        poly=[(x+q[0]*scale,y+q[1]*scale) for q in corners]
        d.line(poly+[poly[0]],fill='#ffdc44',width=3)
        for label,point in landmarks:
            q=G@np.array([*point,1]);px=x+q[0]*scale;py=y+q[1]*scale
            d.ellipse((px-7,py-7,px+7,py+7),outline='#00dcff',width=3)
            d.text((px-34,py-32),label,font=small,fill='#00dcff',stroke_width=1,stroke_fill='black')
    d.text((20,940),'Cyan L1-L3: foil-edge location features, not dislocation labels.',font=font,fill='black')
    d.text((20,978),'Yellow boxes: measured footprints of the 45k images. No defect-type annotations.',font=font,fill='black')
    canvas.save(OUT/'Y5_common_field_localization.png')
    # A larger candidate region centered on a distinctive elongated contrast
    # feature visible in 0031 and 0034. Nature of the feature is unassigned.
    center31=np.array([1040,835,1])
    center28=np.linalg.inv(high[1])@center31
    centers=[M@center28 for M in high]
    canvas=Image.new('RGB',(2490,790),'white');d=ImageDraw.Draw(canvas)
    d.text((20,18),'Y-5.00 | local correspondence check | original-intensity crops',font=font,fill='black')
    d.text((20,53),'Nominal reflection families from conventional alpha-Ti indexing; no Burgers-vector labels.',font=small,fill='black')
    roi=[]
    for j,(aid,im,c) in enumerate(zip(ids,frames,centers)):
        cx,cy=np.rint(c[:2]).astype(int);box=(cx-400,cy-260,cx+400,cy+260)
        x=20+j*825;y=130
        row=['(0002)','{10-10}','{10-11}'][j]
        d.text((x,96),f'{aid} | {row} | center ({cx},{cy})',font=font,fill='black')
        local=im.crop(box)
        local.save(OUT/f'{aid}_local_feature_original_crop.png')
        canvas.paste(local,(x,y))
        # New bar calculated from the original 45k pixel calibration.
        bar=round(100/1.0761216647789286)
        d.line((x+22,y+495,x+22+bar,y+495),fill='white',width=6)
        d.text((x+22,y+459),'100 nm',font=small,fill='white',stroke_width=1,stroke_fill='black')
        roi.append(dict(image_id=aid,box=[int(v) for v in box],selection='elongated-feature neighborhood; type unassigned'))
    d.text((20,698),'The elongated feature in 0031 / 0034 is used for local review, not assigned as a dislocation.',font=font,fill='black')
    d.text((20,738),'The 0028 position is predicted from edge landmarks; absence of a matching trace is not yet an extinction assignment.',font=small,fill='black')
    canvas.save(OUT/'Y5_local_feature_review.png')
    (OUT/'local_feature_roi_coordinates.json').write_text(json.dumps(roi,indent=2)+'\n')


if __name__=='__main__':main()
