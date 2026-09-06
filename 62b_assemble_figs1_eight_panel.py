#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
62b_assemble_figs1_eight_panel.py
Assemble the single-panel PNGs produced by scripts 62 / 38 into Fig. S1 (8 panels, a-h).

Panel composition (finalized 2026-09-05; after Fig. 2 was trimmed to 7 panels,
the batch-correction plots were moved down to SI):
  a  PCA before / Harmony after (colored by dataset, dual combined plot)
        <- path62_batch_PCA_vs_Harmony.png  (rasterized from FigS1_batch_PCA_vs_Harmony.pdf)
  b  blood QC (Genes / UMIs / Mito%, by dataset)   <- path62_blood_qc.png
  c  lung UMAP by cell_type                        <- path62_lung_celltype_umap.png
  d  lung UMAP by cell_category                    <- path62_lung_cell_category_umap.png
  e  lung UMAP by disease                          <- path62_lung_disease_umap.png
  f  lung marker dot plot                          <- path62_lung_dotplot.png
  g  lung QC (Genes / UMIs / Mito%, by disease)    <- path62_lung_qc.png

Layout: each row independently spans the full canvas width; rows scale at equal
height/width ratio; a single panel taller than HMAX is capped by height (to avoid
the dotplot being stretched too large), and the row is centered.
Output: FigS1_blood_lung_QC.{png,pdf}  overwrites Figure/ and server_version/
"""
import os
import sys
from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)

PANELS = [
    ("a", "path62_PCA_before.png"),          # PCA before batch correction (old Fig. 2a content, exported by user)
    ("b", "path62_UMAP_by_dataset.png"),     # cross-dataset mixed UMAP (old Fig. 2d content, exported by user)
    ("c", "path62_blood_qc.png"),
    ("d", "path62_lung_celltype_umap.png"),
    ("e", "path62_lung_cell_category_umap.png"),
    ("f", "path62_lung_disease_umap.png"),
    ("g", "path62_lung_dotplot.png"),
    ("h", "path62_lung_qc.png"),
]
ROWS = [[0, 1, 2], [3, 4, 5], [6, 7]]   # panel indices per row

W_IN    = 17.7      # canvas width (inches)
DPI     = 300
GAP_X   = 55
GAP_Y   = 70
MARGIN  = 45
LABEL_H = 95        # height reserved for panel letter
HMAX    = 1400      # max single-panel height (px)
BG      = (255, 255, 255)

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


def flatten_rgba(im):
    """Composite the alpha channel onto a white background to avoid transparent areas becoming black blocks."""
    if im.mode in ("RGBA", "LA", "P"):
        im = im.convert("RGBA")
        bgw = Image.new("RGBA", im.size, (255, 255, 255, 255))
        im = Image.alpha_composite(bgw, im)
    return im.convert("RGB")


def main():
    imgs = []
    for letter, fname in PANELS:
        f = os.path.join(HERE, fname)
        if not os.path.exists(f):
            sys.exit("[ERROR] missing panel file: %s" % f)
        im = flatten_rgba(Image.open(f))
        imgs.append((letter, im))
        print("  loaded %-38s %dx%d (aspect %.2f)"
              % (fname, im.width, im.height, im.width / float(im.height)))

    W = int(round(W_IN * DPI))

    # ---- compute panel sizes row by row ----
    row_boxes = []          # [ (row_h, [(letter, w, h), ...]) ]
    for row in ROWS:
        n = len(row)
        cw = (W - MARGIN * 2 - GAP_X * (n - 1)) / float(n)
        items, hh = [], []
        for idx in row:
            letter, im = imgs[idx]
            ar = im.width / float(im.height)
            w = int(round(cw))
            h = int(round(cw / ar))
            if h > HMAX:                       # cap by height, shrink width proportionally
                h = HMAX
                w = int(round(HMAX * ar))
            items.append((letter, w, h))
            hh.append(h)
        row_boxes.append((max(hh), items))

    H = MARGIN * 2 + sum(rh + LABEL_H for rh, _ in row_boxes) + GAP_Y * (len(row_boxes) - 1)

    canvas = Image.new("RGB", (W, H), BG)
    dr = ImageDraw.Draw(canvas)
    font = load_font(78)

    y = MARGIN
    for rh, items in row_boxes:
        row_w = sum(w for _, w, _ in items) + GAP_X * (len(items) - 1)
        x = MARGIN + max(0, (W - MARGIN * 2 - row_w) // 2)   # center within row
        for letter, w, h in items:
            idx = [i for i, (L, _) in enumerate(imgs) if L == letter][0]
            im = imgs[idx][1]
            dr.text((x + 4, y), letter, fill=(0, 0, 0), font=font)
            canvas.paste(im.resize((w, h), Image.LANCZOS), (x, y + LABEL_H))
            print("  placed %s -> (%d, %d) size %dx%d" % (letter, x, y + LABEL_H, w, h))
            x += w + GAP_X
        y += LABEL_H + rh + GAP_Y

    out_png = os.path.join(HERE, "FigS1_blood_lung_QC.png")
    out_pdf = os.path.join(HERE, "FigS1_blood_lung_QC.pdf")
    canvas.save(out_png, dpi=(DPI, DPI))
    canvas.save(out_pdf, "PDF", resolution=DPI, save_all=True)
    print("[ok] %s  (%dx%d px = %.2f x %.2f in @%ddpi)"
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
