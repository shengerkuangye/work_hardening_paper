#!/usr/bin/env python3
"""Read-only geometric correspondence audit, not a dislocation classifier.

Dependencies: PYTHONPATH=.codex_tmp/tem_analysis_deps python3 this_script.py
Feature preprocessing is for locating points only. Published/contact images
retain original intensities. Transforms are diagnostics, never new data.
"""
from pathlib import Path
import json
import math
import numpy as np
import cv2
from PIL import Image,ImageDraw,ImageFont

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'data/9_DF/DF'
OUT=ROOT/'results/tem_field_correspondence_2026-10-01'
CACHE=ROOT/'.codex_tmp/tem_correspondence_2026-10-01'
FONT='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
CHANNELS=['BF-S','DF-S','HAADF']
PAIRS=[('0026','0028','within_condition'),('0029','0031','within_condition'),('0032','0034','within_condition'),
       ('0009','0011','within_condition'),('0012','0014','within_condition'),('0015','0017','within_condition'),
       ('0026','0029','cross_condition'),('0029','0032','cross_condition'),('0026','0032','cross_condition'),
       ('0009','0012','cross_condition'),('0012','0015','cross_condition'),('0009','0015','cross_condition')]


def path(aid,ch='BF-S'):
    return next(SRC.glob(aid+'*STEM '+ch+'.tif'))


def features(aid,ch):
    cache=CACHE/f'{aid}_{ch}_features.npz'
    if cache.exists():
        z=np.load(cache);return z['pts'],z['desc']
    with Image.open(path(aid,ch)) as im:arr=np.array(im.convert('L'))
    arr=cv2.resize(arr,(1536,1536),interpolation=cv2.INTER_AREA)
    mask=np.ones(arr.shape,dtype=np.uint8)*255
    mask[1390:,0:550]=0  # exclude scale text/bar, not material contrast
    sift=cv2.SIFT_create(nfeatures=10000,contrastThreshold=.012,edgeThreshold=12)
    key,desc=sift.detectAndCompute(arr,mask)
    pts=np.array([k.pt for k in key],dtype=np.float32)*2048/1536
    np.savez_compressed(cache,pts=pts,desc=desc)
    return pts,desc


def pair(a,b,kind):
    candidate=[]
    for ch in CHANNELS:
        p1,d1=features(a,ch);p2,d2=features(b,ch)
        bf=cv2.BFMatcher()
        m12=bf.knnMatch(d1,d2,k=2);m21=bf.knnMatch(d2,d1,k=2)
        reverse={m.queryIdx:m.trainIdx for m,n in m21 if m.distance<.78*n.distance}
        good=[m for m,n in m12 if m.distance<.75*n.distance and reverse.get(m.trainIdx)==m.queryIdx]
        for m in good:
            candidate.append(dict(channel=ch,p1=p1[m.queryIdx].tolist(),p2=p2[m.trainIdx].tolist(),distance=float(m.distance)))
    result=dict(source=a,target=b,kind=kind,candidate_count=len(candidate),candidates=candidate)
    if len(candidate)<4:return result
    p1=np.float32([m['p1'] for m in candidate]);p2=np.float32([m['p2'] for m in candidate])
    if kind=='within_condition':
        # Nominal magnification change, roughly centered, allowed generous shift.
        scale=2.1703696219522647/1.0761216647789286
        offsets=p2-scale*p1
        plausible=(np.linalg.norm(offsets+1024,axis=1)<1200)
    else:
        # Same 22.5k field should be nearby; excludes unrelated global matches.
        plausible=np.linalg.norm(p2-p1,axis=1)<400
    options=np.where(plausible)[0]
    if len(options)<4:return result
    M,mask=cv2.estimateAffinePartial2D(p1[options],p2[options],method=cv2.RANSAC,
                                    ransacReprojThreshold=9,maxIters=20000,confidence=.999,refineIters=30)
    if M is None:return result
    idx=options[mask.ravel().astype(bool)]
    # Collapse duplicate detections across channels and multiple SIFT scales.
    unique=[]
    for j in idx:
        if not any(np.linalg.norm(p1[j]-p1[k])<18 for k in unique):unique.append(int(j))
    scale=math.hypot(M[0,0],M[1,0]);rotation=math.degrees(math.atan2(M[1,0],M[0,0]))
    expected=2.1703696219522647/1.0761216647789286 if kind=='within_condition' else 1
    geometry_ok=abs(scale/expected-1)<.12 and abs(rotation)<12
    residual=np.linalg.norm(p1[idx]@M[:,:2].T+M[:,2]-p2[idx],axis=1)
    result.update(matrix=M.tolist(),scale=scale,rotation_deg=rotation,
                  inlier_indices=idx.tolist(),unique_inlier_indices=unique,
                  unique_inliers=len(unique),rms_target_pixels=float(np.sqrt(np.mean(residual**2))),
                  geometric_screen_pass=bool(geometry_ok and len(unique)>=6))
    return result


def contact(result):
    if not result.get('geometric_screen_pass'):return
    w=1000;pad=25
    canvas=Image.new('RGB',(2*w+3*pad,w+110),'white');d=ImageDraw.Draw(canvas)
    f=ImageFont.truetype(FONT,24)
    for j,aid in enumerate([result['source'],result['target']]):
        with Image.open(path(aid)) as im:frame=im.convert('RGB').resize((w,w),Image.Resampling.LANCZOS)
        canvas.paste(frame,(pad+j*(w+pad),85))
    d.text((pad,12),f"{result['source']} -> {result['target']} | geometric matches only | {result['unique_inliers']} unique points",font=f,fill='black')
    d.text((pad,43),'Numbers locate matching features; they are not dislocation-type labels.',font=f,fill='black')
    # A sparse, spatially distributed subset is readable; counts are recorded
    # in JSON and are not additional independent experiments.
    selected=[]
    for k in result['unique_inlier_indices']:
        point=np.array(result['candidates'][k]['p1'])
        if not any(np.linalg.norm(point-np.array(result['candidates'][j]['p1']))<260 for j in selected):
            selected.append(k)
        if len(selected)==12:break
    for n,k in enumerate(selected):
        m=result['candidates'][k]
        for j,key in enumerate(['p1','p2']):
            x,y=m[key];x=x*w/2048+pad+j*(w+pad);y=y*w/2048+85
            d.ellipse((x-5,y-5,x+5,y+5),outline='#00cfff',width=2)
            d.text((x+5,y-24),str(n+1),font=f,fill='#00cfff',stroke_width=1,stroke_fill='black')
    canvas.save(OUT/f"matches_{result['source']}_{result['target']}.png")


if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=True);CACHE.mkdir(parents=True,exist_ok=True)
    results=[]
    for a,b,k in PAIRS:
        r=pair(a,b,k);results.append(r);contact(r)
        (OUT/'feature_correspondence_audit.json').write_text(json.dumps(results,indent=2)+'\n')
        print({key:val for key,val in r.items() if key not in ['candidates','inlier_indices','unique_inlier_indices']},flush=True)
