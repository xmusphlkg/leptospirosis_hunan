"""Cache real imagery and administrative reference locations for Figure 1A.

DataV coordinates are GCJ-02; convert them to WGS84 to match Esri imagery.
Markers represent administrative centres, not unrecorded field sampling GPS.
"""
from pathlib import Path
import csv
import gzip
from PIL import Image
import hashlib
import json
import math
import requests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'data/figure1_map_cache'
OUT.mkdir(parents=True, exist_ok=True)

def fetch(url, path, params=None):
    if not path.exists():
        r = requests.get(url, params=params, timeout=90)
        r.raise_for_status()
        path.write_bytes(r.content)
    return path

def wgs(lon, lat):
    x, y = lon - 105, lat - 35
    a = -100 + 2*x + 3*y + .2*y*y + .1*x*y + .2*math.sqrt(abs(x))
    b = 300 + x + 2*y + .1*x*x + .1*x*y + .1*math.sqrt(abs(x))
    a += (20*math.sin(6*x*math.pi)+20*math.sin(2*x*math.pi))*2/3
    b += (20*math.sin(6*x*math.pi)+20*math.sin(2*x*math.pi))*2/3
    a += (20*math.sin(y*math.pi)+40*math.sin(y/3*math.pi))*2/3
    a += (160*math.sin(y/12*math.pi)+320*math.sin(y*math.pi/30))*2/3
    b += (20*math.sin(x*math.pi)+40*math.sin(x/3*math.pi))*2/3
    b += (150*math.sin(x/12*math.pi)+300*math.sin(x/30*math.pi))*2/3
    rad = lat/180*math.pi
    magic = 1-.00669342162296594323*math.sin(rad)**2
    a = a*180/((6378245*(1-.00669342162296594323))/(magic*math.sqrt(magic))*math.pi)
    b = b*180/(6378245/math.sqrt(magic)*math.cos(rad)*math.pi)
    return lon-b, lat-a

base = 'https://geo.datav.aliyun.com/areas_v3/bound/'
bound = json.loads(fetch(base+'430000.json', OUT/'hunan_gcj02.geojson').read_text())
with (OUT/'boundary.csv').open('w') as f:
    writer = csv.writer(f); writer.writerow(['group','lon','lat'])
    for i, polygon in enumerate(bound['features'][0]['geometry']['coordinates']):
        for j, ring in enumerate(polygon):
            for xy in ring: writer.writerow([f'{i}_{j}', *wgs(*xy)])
# Historical labels are retained in the manuscript (Ningxiang is now a city).
sites = {'Ningxiang':430182,'Liuyang':430181,'Xiangtan':430321,
         'Liling':430281,'Junshan':430611,'Yuanjiang':430981,
         'Shuangfeng':431321,'Lianyuan':431382}
with (OUT/'sites.csv').open('w') as f:
    writer = csv.writer(f); writer.writerow(['site','adcode','lon','lat','location_type'])
    for name, code in sites.items():
        obj = json.loads(fetch(base+f'{code}.json',OUT/f'{code}.geojson').read_text())
        writer.writerow([name,code,*wgs(*obj['features'][0]['properties']['center']),'administrative centre; not sampling GPS'])
service = 'https://services.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer'
meta = json.loads(fetch(service, OUT/'imagery_service.json', {'f':'json'}).read_text())
export = json.loads(fetch(service+'/export',OUT/'imagery_export_pretty.json',
    dict(bbox='108,24,115,31',bboxSR=4326,imageSR=4326,size='1800,1800',format='png',f='json')).read_text())
fetch(export['href'],OUT/'satellite_pretty.png')
(OUT/'satellite.png').write_bytes((OUT/'satellite_pretty.png').read_bytes())
with gzip.open(OUT/'satellite.rgb.gz', 'wb') as f:
    f.write(Image.open(OUT/'satellite_pretty.png').convert('RGB').tobytes())
(OUT/'provenance.json').write_text(json.dumps(dict(
    imagery_service=service, attribution=meta.get('copyrightText'),
    boundary_source=base+'430000.json', location_source=base+'{adcode}.json',
    coordinate_note='DataV GCJ-02 converted approximately to WGS84; administrative centres only.',
    retrieved='2026-10-22',extent=export['extent'],
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.name!='provenance.json'}),indent=2))
print(meta.get('copyrightText'))
