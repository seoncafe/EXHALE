import sys
from PIL import Image
import numpy as np
im=np.array(Image.open(sys.argv[1]).convert('RGB')).astype(int)
dark=(im.sum(axis=2)<250)
H,W=dark.shape
# vertical frame lines: columns with long continuous dark runs
def runs(v):
    out=[];s=None
    for i,x in enumerate(v):
        if x and s is None: s=i
        if (not x) and s is not None:
            out.append((s,i-1)); s=None
    if s is not None: out.append((s,len(v)-1))
    return out
hl=[]
for y in range(H):
    for a,b in runs(dark[y]):
        if b-a>300: hl.append((y,a,b))
vl=[]
for x in range(W):
    for a,b in runs(dark[:,x]):
        if b-a>300: vl.append((x,a,b))
print('H',hl); print('V',vl)
