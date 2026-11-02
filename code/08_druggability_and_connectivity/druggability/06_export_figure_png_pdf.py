#!/usr/bin/env python
"""Export the radial drug tree from SVG to a 300-dpi PNG and a vector PDF.

    python src/06_export_figure_png_pdf.py

The SVG itself is written by src/05_render_radial_tree.mjs (d3 + tree_radial.js
under jsdom); this script only converts it, so the three figure files always come
from one layout.

PNG  rasterised with resvg (pip install resvg-py). The SVG is 4000 user units
     wide; at 300 dpi that is a 13.3-inch figure, so the raster is written at
     4000 px and tagged 300 dpi. cairosvg is not usable on this machine -- no
     libcairo DLL is present.
PDF  vector, via svglib + reportlab (no rasterisation, text stays selectable).
     Arial is mapped to Helvetica by reportlab; glyph widths shift very slightly
     against the PNG, which is why the SVG layout is not re-derived here.

Reads   figure/CRC_scTWAS_CellType_Gene_Drug_Indication.svg
Writes  figure/CRC_scTWAS_CellType_Gene_Drug_Indication.png  (300 dpi)
        figure/CRC_scTWAS_CellType_Gene_Drug_Indication.pdf  (vector)
"""
import io
from pathlib import Path

import resvg_py
from PIL import Image
from reportlab.graphics import renderPDF
from svglib.svglib import svg2rlg

FIG = Path(r"W:\VUMC_projects\scTWAS\Drugs\figure")
SVG = FIG / "CRC_scTWAS_CellType_Gene_Drug_Indication.svg"
PNG = SVG.with_suffix(".png")
PDF = SVG.with_suffix(".pdf")

DPI = 300
# Rasterise 1:1 with the SVG's own user units, so the 300-dpi tag is honest and the
# PNG never up- or down-samples the layout (the canvas size changes when rings are
# added or dropped).
import re

PX_WIDTH = round(float(re.search(r'\bwidth="([0-9.]+)"', SVG.read_text(encoding="utf-8")).group(1)))

# --- PNG ---------------------------------------------------------------------
raster = resvg_py.svg_to_bytes(
    svg_path=str(SVG),
    background="#ffffff",
    width=PX_WIDTH,
)
img = Image.open(io.BytesIO(bytes(raster)))
img.save(PNG, dpi=(DPI, DPI))

# --- PDF ---------------------------------------------------------------------
# d3.create("svg") leaves style="max-width: 100%; height: auto" on the root, which
# keeps the figure responsive in a browser but makes svglib try to read a length
# from "auto". Drop it on a scratch copy; the file on disk keeps the attribute.
import tempfile
import xml.etree.ElementTree as ET

ET.register_namespace("", "http://www.w3.org/2000/svg")
tree = ET.parse(SVG)
tree.getroot().attrib.pop("style", None)
with tempfile.TemporaryDirectory() as tmp:
    flat = Path(tmp) / SVG.name
    tree.write(flat, encoding="utf-8", xml_declaration=True)
    drawing = svg2rlg(str(flat))

renderPDF.drawToFile(drawing, str(PDF))

print(f"PNG {img.size[0]}x{img.size[1]} px @ {DPI} dpi "
      f"({img.size[0] / DPI:.1f} x {img.size[1] / DPI:.1f} in), "
      f"{PNG.stat().st_size / 1e6:.1f} MB")
print(f"PDF {drawing.width:.0f}x{drawing.height:.0f} pt, "
      f"{PDF.stat().st_size / 1e6:.1f} MB")
