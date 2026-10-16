"""Print tick-mark pixel positions next to a frame edge.
usage: ticks.py img side pos lo hi   (side: bottom/top -> pos is row; left/right -> pos is column)"""
import sys
import numpy as np
from PIL import Image
im = np.asarray(Image.open(sys.argv[1]).convert('RGB')).astype(int)
side, pos, lo, hi = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
# dark or colored (non-white) pixels
nonwhite = (im.min(axis=2) < 150)
def groups(idx):
    out=[]; 
    for i in idx:
        if out and i-out[-1][-1] <= 1: out[-1].append(i)
        else: out.append([i])
    return [ (g[0]+g[-1])/2 for g in out]
for off in (range(2,5), range(-4,-1)):
    if side in ('bottom','top'):
        band = nonwhite[[pos+o for o in off], lo:hi]
        hit = np.where(band.all(axis=0))[0] + lo
    else:
        band = nonwhite[lo:hi, [pos+o for o in off]]
        hit = np.where(band.all(axis=1))[0] + lo
    print('offset', list(off)[0], groups(hit))
