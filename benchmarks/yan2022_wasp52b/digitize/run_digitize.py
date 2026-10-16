import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from digitize_yan2022 import *
import numpy as np
from PIL import Image, ImageDraw

# directory holding the 300 dpi page renders hi15-15.png, hi18-18.png, hi29-29.png:
#   pdftoppm -f N -l N -r 300 -png 'references/Yan 2022.pdf' hiN
D = os.environ.get('YAN_RENDER_DIR', os.getcwd())
specs = {
 # name: (png, frame, ycal, ylog, ntrack)
 'fig9a': ('hi15-15.png', (377.5, 1267.5, 271, 934), (929.5, 0.0, 0.0576), False, 1),
 'fig9b': ('hi15-15.png', (1404.5, 2291, 274, 935), (904.5, 0.0, 8.4), False, 1),
 'fig9c': ('hi15-15.png', (349, 1207.5, 1184.5, 1824.5), (1792.0, 1e3, 71.64), True, 2),
 'fig23a': ('hi29-29.png', (377.5, 1267.5, 271, 934), (922.5, 0.0, 0.0583), False, 1),
 'fig23b': ('hi29-29.png', (1404.5, 2291, 274, 935), (904.3, 0.0, 8.325), False, 1),
 'fig23c': ('hi29-29.png', (337, 1207, 1185, 1833.5), (1833.5, 10.0, 68.5), True, 1),
 'fig23d': ('hi29-29.png', (1375.5, 2231.5, 1197.5, 1835), (1835.0, 1e-3, 58.0), True, 2),
 'fig12c': ('hi18-18.png', (349, 1207.5, 1184.5, 1824.5), (1792.0, 1e3, 72.79), True, 2),
}
res = {}
ims = {}
for name, (png, frame, ycal, ylog, nt) in specs.items():
    if png not in ims:
        ims[png] = load(os.path.join(D, png))
    cols = ['red'] if name == 'fig12c' else ['0.5', '1', '2']
    for cn in cols:
        ex = [(1000, 1210, 1185, 1450)] if name == 'fig23c' else []
        tr = extract(ims[png], frame, COL[cn], ntrack=nt, exclude=ex)
        for k, t in enumerate(tr):
            res[(name, cn, k)] = (t, to_data(t, frame, ycal, ylog))
# overlay check images
for png in ims:
    img = Image.open(os.path.join(D, png)).convert('RGB'); dr = ImageDraw.Draw(img)
    for (name, cn, k), (t, _) in res.items():
        if specs[name][0] != png: continue
        for x, y in t[::4]:
            dr.ellipse((x-2, y-2, x+2, y+2), fill=(0, 0, 0) if k == 0 else (120, 60, 0))
    img.save(os.path.join(D, 'check_' + png))
import pickle
pickle.dump({k: v[1] for k, v in res.items()}, open(os.path.join(D, 'digitized.pkl'), 'wb'))
for k, v in res.items():
    d = v[1]
    if d:
        print(k, len(d), 'r %.3f-%.3f' % (d[0][0], d[-1][0]))
    else:
        print(k, 0)
