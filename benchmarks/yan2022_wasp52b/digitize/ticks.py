import sys
from PIL import Image
import numpy as np
f=sys.argv[1]; x0,x1,y0,y1=[float(v) for v in sys.argv[2:6]]
im=np.array(Image.open(f).convert('RGB')).astype(int)
dark=(im.sum(axis=2)<300)
x0i,x1i,y0i,y1i=int(round(x0)),int(round(x1)),int(round(y0)),int(round(y1))
# left ticks: outside the frame, columns x0-10..x0-4
def groups(idx):
    g=[];cur=[]
    for i in idx:
        if cur and i-cur[-1]>1: g.append(cur);cur=[]
        cur.append(i)
    if cur: g.append(cur)
    return [ (np.mean(c),len(c)) for c in g]
col=dark[y0i-2:y1i+3, x0i+2:x0i+6].all(axis=1)
print('left ticks (y):',[ (round(m+y0i-2,1),n) for m,n in groups(np.where(col)[0])])
row=dark[y1i-6:y1i-2, x0i-2:x1i+3].all(axis=0)
print('bottom ticks (x):',[ (round(m+x0i-2,1),n) for m,n in groups(np.where(row)[0])])
