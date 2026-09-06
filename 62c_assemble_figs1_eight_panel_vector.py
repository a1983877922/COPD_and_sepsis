#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
62c_assemble_figs1_eight_panel_vector.py

Assemble the vector PDFs produced by scripts 62/38 plus the S1A/S1B panels exported by
the user from Illustrator into a **vector** layout of Fig. S1 (8 panels a-h) — every
element remains selectable/editable in Illustrator.

Key: pymupdf.show_pdf_page() pastes the source PDF page onto the target rectangle as a
Form XObject, preserving all vector paths/text as-is; panel letters are drawn with
pymupdf.insert_text() using the Helvetica-Bold vector font.

Panels (finalized 2026-09-05 18:30):
  a  PCA before        <- Figure/FigS1A.…pdf (exported by user from Illustrator)
  b  UMAP by dataset   <- Figure/FigS1B.…pdf (exported by user from Illustrator)
  c  blood QC          <- server_version/path62_blood_qc.pdf
  d  lung celltype     <- server_version/path62_lung_celltype_umap.pdf
  e  lung category     <- server_version/path62_lung_cell_category_umap.pdf
  f  lung disease      <- server_version/path62_lung_disease_umap.pdf
  g  lung dotplot      <- server_version/path62_lung_dotplot.pdf
  h  lung QC           <- server_version/path62_lung_qc.pdf

Output: FigS1_blood_lung_QC.pdf  overwrites server_version/ and Figure/
      (also writes a _vector.pdf copy for provenance)
"""
import os
import sys
import shutil
import pymupdf

HERE   = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
FIGDIR = os.path.join(PARENT, "Figure")

PANELS = [
    ("a", os.path.join(FIGDIR,  "FigS1A.PCA of the integrated atlas before batch correction, coloured by dataset of origin; cells separate strongly by dataset..pdf")),
    ("b", os.path.join(FIGDIR,  "FigS1B.UMAP coloured by dataset of origin, showing cross-dataset mixing..pdf")),
    ("c", os.path.join(HERE,    "path62_blood_qc.pdf")),
    ("d", os.path.join(HERE,    "path62_lung_celltype_umap.pdf")),
    ("e", os.path.join(HERE,    "path62_lung_cell_category_umap.pdf")),
    ("f", os.path.join(HERE,    "path62_lung_disease_umap.pdf")),
    ("g", os.path.join(HERE,    "path62_lung_dotplot.pdf")),
    ("h", os.path.join(HERE,    "path62_lung_qc.pdf")),
]
ROWS = [(0, 1, 2), (3, 4, 5), (6, 7)]          # panel indices per row

COLS, ROWS_N  = 3, 3                           # 3x3 layout
PANEL_W       = 3.4                            # single-panel width (in)
PANEL_H       = 3.0                            # single-panel height (in) — scaled by HMAX
LABEL_H       = 0.35                           # letter row height (in)
GAP_X, GAP_Y  = 0.18, 0.20
MARGIN        = 0.30
FONT_NAME     = "Helvetica-Bold"
FONT_SIZE     = 18                             # pt


def panel_size(src_pdf):
    """Read the physical size of page 0 of the source PDF, return (w_in, h_in)."""
    d = pymupdf.open(src_pdf)
    r = d[0].rect
    d.close()
    return r.width / 72.0, r.height / 72.0


def fit_rect(src_pdf, target_w_in, target_h_in):
    """Fit the source panel proportionally into the target box; return the fitted
    rectangle (x0,y0,x1,y1) in pt — centered."""
    sw, sh = panel_size(src_pdf)
    ar_s = sw / sh
    ar_t = target_w_in / target_h_in
    if ar_s >= ar_t:                            # source wider, scale by width
        w = target_w_in
        h = w / ar_s
    else:                                       # source taller, scale by height
        h = target_h_in
        w = h * ar_s
    x0 = (target_w_in - w) / 2
    y0 = (target_h_in - h) / 2
    return x0 * 72, y0 * 72, (x0 + w) * 72, (y0 + h) * 72


def main():
    # ---- 1. verify all sources exist ----
    for letter, src in PANELS:
        if not os.path.exists(src):
            sys.exit("[ERROR] missing panel: %s" % src)
        w, h = panel_size(src)
        print("  %-3s  %5.2fx%5.2f in  <-  %s" % (letter, w, h, os.path.basename(src)))

    # ---- 2. create empty canvas ----
    page_w_in = MARGIN * 2 + COLS * PANEL_W + (COLS - 1) * GAP_X
    page_h_in = MARGIN * 2 + ROWS_N * (PANEL_H + LABEL_H) + (ROWS_N - 1) * GAP_Y
    page_w_pt = page_w_in * 72
    page_h_pt = page_h_in * 72
    print("\n[canvas] %.2f x %.2f in (%.0f x %.0f pt)" %
          (page_w_in, page_h_in, page_w_pt, page_h_pt))

    out = pymupdf.open()
    page = out.new_page(width=page_w_pt, height=page_h_pt)

    # ---- 3. paste vector content panel by panel ----
    for ridx, row in enumerate(ROWS):
        for cidx, pidx in enumerate(row):
            letter, src = PANELS[pidx]

            # target box (in -> pt); bottom-left is the origin, pymupdf origin is top-left
            # pymupdf coords: y0 runs downward from the top of the page, so y_top counts from the top
            y_top = MARGIN + ridx * (PANEL_H + LABEL_H + GAP_Y)
            y_bot = y_top + PANEL_H + LABEL_H
            x_lft = MARGIN + cidx * (PANEL_W + GAP_X)
            x_rgt = x_lft + PANEL_W

            # letter label position: 0.05 in padding below the panel box
            label_x = x_lft * 72 + 2
            label_y = y_top * 72 + 2
            page.insert_text(
                (label_x, label_y + 18),
                letter,
                fontname=FONT_NAME,
                fontsize=FONT_SIZE,
                color=(0, 0, 0),
            )

            # vector paste (Form XObject, editable in Illustrator)
            src_doc = pymupdf.open(src)
            # fit rectangle: (x0,y0,x1,y1) in pt, scaled proportionally to the source page physical size
            fit_rect_pt = fit_rect(src, PANEL_W, PANEL_H)
            # flip y into pymupdf coordinate system (top-down)
            target_rect = pymupdf.Rect(
                x_lft * 72 + fit_rect_pt[0],
                y_top * 72 + LABEL_H * 72 + fit_rect_pt[1],
                x_lft * 72 + fit_rect_pt[2],
                y_top * 72 + LABEL_H * 72 + fit_rect_pt[3],
            )
            page.show_pdf_page(target_rect, src_doc, 0, keep_proportion=False)
            src_doc.close()
            print("  placed %s @ (%.2f, %.2f) %.2fx%.2f in -> fit %s" %
                  (letter, x_lft, y_top + LABEL_H, PANEL_W, PANEL_H, fit_rect_pt))

    # ---- 4. output ----
    out_pdf = os.path.join(HERE, "FigS1_blood_lung_QC.pdf")
    out.save(out_pdf, deflate=True, garbage=4)
    out.close()
    print("\n[ok] %s  (vector, Illustrator-editable)" % out_pdf)

    # sync to Figure/
    if os.path.isdir(FIGDIR):
        shutil.copy2(out_pdf, os.path.join(FIGDIR, "FigS1_blood_lung_QC.pdf"))
        print("[ok] synced -> %s" % os.path.join(FIGDIR, "FigS1_blood_lung_QC.pdf"))


if __name__ == "__main__":
    main()
