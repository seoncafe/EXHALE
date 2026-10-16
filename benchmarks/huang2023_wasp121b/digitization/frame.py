"""Locate the plot frame and tick marks in a cropped 300 dpi figure image."""
import sys
import numpy as np
from PIL import Image
im = np.asarray(Image.open(sys.argv[1]).convert('RGB')).astype(int)
dark = (im.sum(axis=2) < 200)
H, W = dark.shape
rows = dark.sum(axis=1); cols = dark.sum(axis=0)
print('size', W, H)
print('long rows', [(i, rows[i]) for i in range(H) if rows[i] > 0.5*W])
print('long cols', [(j, cols[j]) for j in range(W) if cols[j] > 0.4*H])
