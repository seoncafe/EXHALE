"""Overlay traced pixel tracks on the figure image for a visual check."""
import sys
import numpy as np
from PIL import Image, ImageDraw
def overlay(img, tracks, out):
    im = Image.open(img).convert('RGB')
    d = ImageDraw.Draw(im)
    colors = [(0, 200, 0), (255, 140, 0), (0, 200, 255), (200, 0, 200), (120, 120, 0), (0, 0, 0)]
    for k, (cols, rows) in enumerate(tracks):
        for c, r in zip(cols, rows):
            d.ellipse((c - 1.5, r - 1.5, c + 1.5, r + 1.5), outline=colors[k % 6])
    im.save(out)
