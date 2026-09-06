#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
62b_assemble_figs1_six_panel.py
Assemble the single-panel PNGs produced by script 62 into Fig. S1 (6 panels, a-f).

Panel composition (per user decision on 2026-09-05):
  a  blood QC (genes/UMIs/Mito, by dataset)        <- path62_blood_qc.png
  b  lung UMAP by cell_type                        <- path62_lung_celltype_umap.png
  c  lung UMAP by cell_category                    <- path62_lung_cell_category_umap.png
  d  lung UMAP by disease                          <- path62_lung_disease_umap.png
  e  lung marker dot plot                          <- path62_lung_dotplot.png
  f  lung QC (genes/UMIs/Mito, by disease)         <- path62_lung_qc.png

(The original blood-side a-d: celltype UMAP / dotplot / group UMAP / dataset UMAP
 were removed because they duplicated Fig. 2c/2g/2e/2d one-to-one.)

Layout: 3 columns x 2 rows, equal-width columns (so all panels share the same font
size), vertically centered within each row.
Output: FigS1_blood_lung_QC.{png,pdf}  overwrites Figure/ and server_version/
"""
import os
import sys
from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)

PANELS = [
    ("a", "path62_blood_qc.png"),
    ("b", "path62_lung_celltype_umap.png"),
    ("c", "path62_lung_cell_category_umap.png"),
    ("d", "path62_lung_disease_umap.png"),
    ("e", "path62_lung_dotplot.png"),
    ("f", "path62_lung_qc.png"),
]

COL_W    = 1700   # target width per column (px)
GAP_X    = 55
GAP_Y    = 70
MARGIN   = 45
LABEL_H  = 95     # height reserved for panel letter
BG       = (255, 255, 255)
DPI      = 300

FONT_CANDIDATES = [
    r"C:\Windows\Fonts\arialbd.ttf",
    r"C:\Windows\Fonts\DejaVuSans-Bold.ttf",
]


def load_font(size):
    for p in FONT_CANDIDATES:
        if os.path.exists(p):
            try:
                return ImageFont.truetype(p, size)
            except Exception:
                pass
    return ImageFont.load_default()


def main():
    imgs, sizes = [], []
    for letter, fname in PANELS:
        f = os.path.join(HERE, fname)
        if not os.path.exists(f):
            sys.exit("[ERROR] missing panel file: %s" % f)
        im = Image.open(f)
        if im.mode in ("RGBA", "LA", "P"):
            im = im.convert("RGBA")
            bgw = Image.new("RGBA", im.size, (255, 255, 255, 255))
            im = Image.alpha_composite(bgw, im)
        im = im.convert("RGB")
        imgs.append((letter, im))
        sizes.append((im.width, im.height))
        print("  loaded %-38s %dx%d" % (fname, im.width, im.height))

    def scaled(w0, h0):
        w = COL_W
        h = int(round(h0 * COL_W / float(w0)))
        return w, h

    sc = [scaled(w0, h0) for (w0, h0) in sizes]
    row_h = [max(sc[0][1], sc[1][1], sc[2][1]),
             max(sc[3][1], sc[4][1], sc[5][1])]

    W = MARGIN * 2 + COL_W * 3 + GAP_X * 2
    H = MARGIN * 2 + LABEL_H * 2 + row_h[0] + GAP_Y + row_h[1]

    canvas = Image.new("RGB", (W, H), BG)
    dr = ImageDraw.Draw(canvas)
    font = load_font(78)

    for i, (letter, im) in enumerate(imgs):
        row, col = divmod(i, 3)
        w, h = sc[i]
        im2 = im.resize((w, h), Image.LANCZOS)
        x = MARGIN + col * (COL_W + GAP_X)
        y = MARGIN + row * (LABEL_H + row_h[row] + GAP_Y)
        dr.text((x + 4, y - 90), letter, fill=(0, 0, 0), font=font)
        y_img = y + (row_h[row] - h) // 2
        canvas.paste(im2, (x, y_img))
        print("  placed %s -> (%d, %d) size %dx%d" % (letter, x, y_img, w, h))

    out_png = os.path.join(HERE, "FigS1_blood_lung_QC.png")
    out_pdf = os.path.join(HERE, "FigS1_blood_lung_QC.pdf")
    canvas.save(out_png, dpi=(DPI, DPI))
    canvas.save(out_pdf, "PDF", resolution=DPI, save_all=True)
    print("[ok] %s  (%dx%d px = %.1f x %.1f in @%ddpi)"
          % (out_png, W, H, W / float(DPI), H / float(DPI), DPI))
    print("[ok] %s" % out_pdf)

    dst = os.path.join(PARENT, "Figure")
    if os.path.isdir(dst):
        import shutil
        for p in (out_png, out_pdf):
            shutil.copy2(p, os.path.join(dst, os.path.basename(p)))
        print("[ok] synced -> %s" % dst)


if __name__ == "__main__":
    main()
