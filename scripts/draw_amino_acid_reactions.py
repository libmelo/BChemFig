"""Generate textbook-ready amino-acid reaction schemes as SVG and PNG."""

import argparse
from html import escape
from io import BytesIO
from math import cos, sin, pi
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent
OUT = Path.cwd() / "amino-acid-reactions"


def find_rdkit_vendor():
    candidates = []
    configured = os.environ.get("RDKIT_VENDOR_PATH")
    if configured:
        candidates.append(Path(configured).expanduser())
    for parent in (ROOT, *ROOT.parents):
        candidates.append(parent / ".vendor" / "rdkit")
    return next((path.resolve() for path in candidates if path.is_dir()), None)


VENDOR = find_rdkit_vendor()
if VENDOR:
    sys.path.insert(0, str(VENDOR))
    if hasattr(os, "add_dll_directory"):
        for folder in (VENDOR / "bin", VENDOR / "rdkit.libs", VENDOR / "numpy.libs"):
            if folder.is_dir():
                os.add_dll_directory(str(folder))

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from PIL import ImageFont


INK = "#111111"
CHEM_FONT = "Times New Roman"
CN_FONT = "SimSun"
CHEM_FONT_FILE = Path(r"C:\Windows\Fonts\times.ttf")
CN_FONT_FILE = Path(r"C:\Windows\Fonts\simsun.ttc")
LINE = 4
_FONT_CACHE = {}


def text_font(path, size):
    key = (str(path), int(size))
    if key not in _FONT_CACHE:
        _FONT_CACHE[key] = ImageFont.truetype(str(path), int(size))
    return _FONT_CACHE[key]


class SVG:
    def __init__(self, width, height, white=True):
        self.width = width
        self.height = height
        self.parts = [
            '<?xml version="1.0" encoding="UTF-8"?>',
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        ]
        if white:
            self.parts.append(f'<rect width="{width}" height="{height}" fill="white"/>')

    def line(self, x1, y1, x2, y2, width=LINE):
        self.parts.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{INK}" stroke-width="{width}" stroke-linecap="round"/>'
        )

    def polygon(self, points, fill=INK):
        pts = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
        self.parts.append(f'<polygon points="{pts}" fill="{fill}"/>')

    def polyline(self, points, width=LINE, close=False):
        pts = list(points)
        if close:
            pts.append(pts[0])
        data = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
        self.parts.append(
            f'<polyline points="{data}" fill="none" stroke="{INK}" '
            f'stroke-width="{width}" stroke-linejoin="round" stroke-linecap="round"/>'
        )

    def text(self, x, y, value, size=54, anchor="middle", family=CHEM_FONT,
             italic=False, bold=False):
        style = "italic" if italic else "normal"
        weight = "700" if bold else "400"
        self.parts.append(
            f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
            f'dominant-baseline="middle" font-family="{family}" font-size="{size}" '
            f'font-style="{style}" font-weight="{weight}" fill="{INK}">{escape(value)}</text>'
        )

    def mixed_text(self, x, y, value, size=40, anchor="middle", bold=False):
        def chinese(ch):
            return "\u3400" <= ch <= "\u9fff" or ch in "（），。：；、“”"

        groups = []
        current = value[0]
        current_cn = chinese(value[0])
        for ch in value[1:]:
            is_cn = chinese(ch)
            if is_cn == current_cn:
                current += ch
            else:
                groups.append((current, current_cn))
                current, current_cn = ch, is_cn
        groups.append((current, current_cn))
        weight = "700" if bold else "400"
        widths = [
            text_font(CN_FONT_FILE if is_cn else CHEM_FONT_FILE, size).getlength(part)
            for part, is_cn in groups
        ]
        total_width = sum(widths)
        if anchor == "middle":
            cursor = x - total_width / 2
        elif anchor == "end":
            cursor = x - total_width
        else:
            cursor = x
        for (part, is_cn), width in zip(groups, widths):
            family = CN_FONT if is_cn else CHEM_FONT
            self.parts.append(
                f'<text x="{cursor:.1f}" y="{y:.1f}" text-anchor="start" '
                f'dominant-baseline="middle" font-family="{family}" font-size="{size}" '
                f'font-weight="{weight}" fill="{INK}">{escape(part)}</text>'
            )
            cursor += width

    def finish(self):
        return "\n".join(self.parts + ["</svg>"])


def arrow(s, x1, x2, y, label_top=None, label_bottom=None):
    s.line(x1, y, x2 - 22, y, 4)
    s.polygon([(x2, y), (x2 - 25, y - 13), (x2 - 25, y + 13)])
    if label_top:
        s.mixed_text((x1 + x2) / 2, y - 51, label_top, 39, bold=True)
    if label_bottom:
        s.text((x1 + x2) / 2, y + 50, label_bottom, 38, bold=True)


def equilibrium_arrow(s, x1, x2, y, top="−H₂O", bottom="+H₂O"):
    s.line(x1, y - 11, x2 - 22, y - 11, 3.5)
    s.polygon([(x2, y - 11), (x2 - 24, y - 23), (x2 - 24, y + 1)])
    s.line(x1 + 22, y + 16, x2, y + 16, 3.5)
    s.polygon([(x1, y + 16), (x1 + 24, y + 4), (x1 + 24, y + 28)])
    s.text((x1 + x2) / 2, y - 61, top, 38, bold=True)
    s.text((x1 + x2) / 2, y + 65, bottom, 38, bold=True)


def plus(s, x, y):
    s.text(x, y, "+", 55)


def caption(s, x, y, value, size=43):
    s.mixed_text(x, y, value, size, bold=True)


def draw_r_prime(s, x, y, size=54):
    """Draw R′ with a consistently positioned superscript prime."""
    s.text(x - 7, y, "R", size)
    s.text(x + 21, y - 18, "′", int(size * 0.58))


def draw_r_prime_x(s, x, y, size=57):
    s.text(x - 34, y, "R", size)
    s.text(x - 7, y - 19, "′", int(size * 0.58))
    s.text(x + 31, y, "X", size)


def draw_r_prime_hn_left(s, end_x, y, size=53):
    """Draw R′HN so the bond at the right terminates at N."""
    s.text(end_x - 96, y, "R", size)
    s.text(end_x - 70, y - 18, "′", int(size * 0.58))
    s.text(end_x - 55, y, "H", size)
    s.text(end_x - 18, y, "N", size)


def double_bond(s, x1, y1, x2, y2, gap=7, width=3.5):
    dx, dy = x2 - x1, y2 - y1
    length = (dx * dx + dy * dy) ** 0.5
    ox, oy = -dy / length * gap / 2, dx / length * gap / 2
    s.line(x1 + ox, y1 + oy, x2 + ox, y2 + oy, width)
    s.line(x1 - ox, y1 - oy, x2 - ox, y2 - oy, width)


def draw_benzene(s, cx, cy, radius=90):
    pts = [(cx + radius * cos(i * pi / 3), cy + radius * sin(i * pi / 3)) for i in range(6)]
    s.polyline(pts, close=True)
    for i, j in ((0, 1), (2, 3), (4, 5)):
        p1, p2 = pts[i], pts[j]
        q1 = (p1[0] * 0.84 + cx * 0.16, p1[1] * 0.84 + cy * 0.16)
        q2 = (p2[0] * 0.84 + cx * 0.16, p2[1] * 0.84 + cy * 0.16)
        s.line(q1[0], q1[1], q2[0], q2[1], 3)
    return pts


def radial_substituent(s, cx, cy, point, label, length=70, size=49, family=CHEM_FONT):
    vx, vy = point[0] - cx, point[1] - cy
    mag = (vx * vx + vy * vy) ** 0.5
    ux, uy = vx / mag, vy / mag
    end = (point[0] + ux * length, point[1] + uy * length)
    s.line(point[0], point[1], end[0], end[1])
    tx, ty = end[0] + ux * 25, end[1] + uy * 25
    anchor = "start" if ux > 0.35 else "end" if ux < -0.35 else "middle"
    s.text(tx, ty, label, size, anchor=anchor, family=family)
    return end


def draw_vertical_connected_label(s, x, y, label, size=53):
    """Place the atom carrying the vertical bond exactly on the bond axis."""
    attached = {
        "COOH": ("C", "OOH"),
        "NH₂": ("N", "H₂"),
        "OH": ("O", "H"),
    }
    if label in attached:
        atom, rest = attached[label]
        s.text(x, y, atom, size)
        s.text(x + 22, y, rest, size, anchor="start")
    elif label == "R":
        s.text(x, y, label, size)
    else:
        s.text(x, y, label, size)


def draw_cross_projection(s, x, y, left, right, top, bottom):
    """Projection with every bond terminating at the chemically attached atom."""
    s.text(x, y, "C", 58)
    if left == "R′HN":
        draw_r_prime_hn_left(s, x - 88, y, 53)
    else:
        s.text(x - 88, y, left, 53, anchor="end")
    s.line(x - 84, y, x - 27, y)
    s.text(x + 88, y, right, 53, anchor="start")
    s.line(x + 27, y, x + 84, y)
    s.line(x, y - 31, x, y - 74)
    draw_vertical_connected_label(s, x, y - 104, top, 53)
    s.line(x, y + 31, x, y + 74)
    draw_vertical_connected_label(s, x, y + 104, bottom, 53)


def draw_horizontal_aa(s, x, y, left="H₂N", top="R", right="COOH", center="CH"):
    draw_cross_projection(s, x, y, left, right, top, "H")


def draw_fischer_aa(s, x, y, left="H₂N", right="H", top="COOH", bottom="R"):
    draw_cross_projection(s, x, y, left, right, top, bottom)


def draw_vertical_substituted_aa(s, x, y, top_group, left="R", right="COOH"):
    draw_cross_projection(s, x, y, left, right, top_group, "H")


def scheme_nitrous(s):
    y = 250
    draw_vertical_substituted_aa(s, 300, y, "NH₂")
    plus(s, 580, y)
    s.text(690, y, "HNO₂", 56)
    arrow(s, 800, 920, y)
    draw_vertical_substituted_aa(s, 1135, y, "OH")
    plus(s, 1420, y)
    s.text(1515, y, "N₂↑", 56)
    plus(s, 1630, y)
    s.text(1735, y, "H₂O", 56)
    caption(s, 300, 440, "α-氨基酸")
    caption(s, 1135, 440, "α-羟基酸")


def scheme_acylation(s):
    y = 270
    draw_fischer_aa(s, 300, y)
    plus(s, 540, y)
    draw_r_prime_x(s, 640, y, 57)
    arrow(s, 750, 870, y)
    draw_fischer_aa(s, 1120, y, left="R′HN")
    plus(s, 1370, y)
    s.text(1460, y, "HX", 56)
    caption(s, 300, 470, "氨基酸")
    caption(s, 1120, 470, "N-酰基氨基酸")


def draw_dnfb_ring(s, cx, cy, product=False):
    pts = draw_benzene(s, cx, cy, 92)
    radial_substituent(s, cx, cy, pts[3], "O₂N", 62, 48)
    radial_substituent(s, cx, cy, pts[5], "NO₂", 64, 48)
    if not product:
        radial_substituent(s, cx, cy, pts[0], "F", 62, 52)
    else:
        nx = pts[0][0] + 84
        s.line(pts[0][0], pts[0][1], nx - 27, pts[0][1])
        s.text(nx, pts[0][1], "N", 50)
        s.line(nx, pts[0][1] - 28, nx, pts[0][1] - 65)
        s.text(nx, pts[0][1] - 91, "H", 45)
        return nx, pts[0][1]
    return None


def scheme_dnfb(s):
    y = 270
    draw_dnfb_ring(s, 300, y, False)
    plus(s, 555, y)
    draw_horizontal_aa(s, 800, y)
    arrow(s, 1075, 1160, y, label_top="弱碱中")
    nx, ny = draw_dnfb_ring(s, 1470, y, True)
    ca = nx + 105
    s.line(nx + 27, ny, ca - 27, ny)
    s.text(ca, ny, "C", 56)
    s.line(ca, ny - 31, ca, ny - 74)
    s.text(ca, ny - 104, "R", 52)
    s.line(ca, ny + 31, ca, ny + 74)
    s.text(ca, ny + 104, "H", 49)
    s.line(ca + 27, ny, ca + 82, ny)
    s.text(ca + 86, ny, "COOH", 53, anchor="start")
    plus(s, 2130, y)
    s.text(2220, y, "HF", 54)
    caption(s, 300, 480, "1-氟-2,4-二硝基苯（DNFB）", 41)
    caption(s, 1640, 480, "DNP-氨基酸（黄色）", 42)


def draw_pitc(s, cx, cy):
    pts = draw_benzene(s, cx, cy, 78)
    s.line(pts[0][0], pts[0][1], pts[0][0] + 52, pts[0][1])
    nx = pts[0][0] + 83
    s.text(nx, cy, "N", 51)
    cx2 = nx + 105
    double_bond(s, nx + 35, cy, cx2 - 35, cy)
    s.text(cx2, cy, "C", 51)
    sx = cx2 + 108
    double_bond(s, cx2 + 35, cy, sx - 35, cy)
    s.text(sx, cy, "S", 51)


def draw_ptc(s, cx, cy):
    pts = draw_benzene(s, cx, cy, 72)
    s.line(pts[0][0], pts[0][1], pts[0][0] + 48, pts[0][1])
    n1x = pts[0][0] + 83
    s.text(n1x, cy, "N", 48)
    s.line(n1x, cy - 27, n1x, cy - 61)
    s.text(n1x, cy - 85, "H", 43)
    cth = n1x + 116
    s.line(n1x + 27, cy, cth - 32, cy)
    s.text(cth, cy, "C", 49)
    s.line(cth, cy + 35, cth, cy + 83)
    s.line(cth + 9, cy + 35, cth + 9, cy + 83)
    s.text(cth + 5, cy + 118, "S", 49)
    n2x = cth + 116
    s.line(cth + 34, cy, n2x - 40, cy)
    s.text(n2x, cy, "N", 48)
    s.line(n2x, cy - 27, n2x, cy - 61)
    s.text(n2x, cy - 85, "H", 43)
    chx = n2x + 126
    s.line(n2x + 27, cy, chx - 27, cy)
    s.text(chx, cy, "C", 50)
    s.line(chx, cy - 29, chx, cy - 72)
    s.text(chx, cy - 100, "R", 47)
    s.line(chx, cy + 29, chx, cy + 72)
    s.text(chx, cy + 100, "H", 44)
    s.line(chx + 27, cy, chx + 82, cy)
    s.text(chx + 86, cy, "COOH", 49, anchor="start")


def draw_pth(s, cx, cy):
    # Regularized 3-phenyl-2-thioxo-5-substituted imidazolidin-4-one ring.
    ring = draw_benzene(s, cx - 255, cy - 60, 67)
    n1 = (cx - 90, cy - 60)
    c_o = (cx + 70, cy - 60)
    c_alpha = (cx + 70, cy + 80)
    n_h = (cx - 10, cy + 165)
    c_s = (cx - 90, cy + 80)

    s.line(ring[0][0], ring[0][1], n1[0] - 31, n1[1])
    s.text(n1[0], n1[1], "N", 48)
    s.line(n1[0] + 31, n1[1], c_o[0] - 29, c_o[1])
    s.text(c_o[0], c_o[1], "C", 47)
    s.line(c_o[0], c_o[1] + 31, c_alpha[0], c_alpha[1] - 31)
    s.text(c_alpha[0], c_alpha[1], "C", 47)
    s.line(c_alpha[0] - 21, c_alpha[1] + 27, n_h[0] + 27, n_h[1] - 24)
    s.text(n_h[0], n_h[1], "N", 46)
    s.line(n_h[0], n_h[1] + 27, n_h[0], n_h[1] + 58)
    s.text(n_h[0], n_h[1] + 82, "H", 41)
    s.line(n_h[0] - 27, n_h[1] - 24, c_s[0] + 21, c_s[1] + 27)
    s.text(c_s[0], c_s[1], "C", 47)
    s.line(c_s[0], c_s[1] - 31, n1[0], n1[1] + 31)

    # Carbonyl and thiocarbonyl are drawn outside the regular five-membered ring.
    double_bond(s, c_o[0], c_o[1] - 32, c_o[0], c_o[1] - 86, 8, 3)
    s.text(c_o[0], c_o[1] - 116, "O", 47)
    double_bond(s, c_s[0], c_s[1] + 32, c_s[0], c_s[1] + 84, 8, 3)
    s.text(c_s[0], c_s[1] + 115, "S", 47)

    # The alpha-carbon is labeled C because H is shown explicitly below it.
    s.line(c_alpha[0] + 32, c_alpha[1], c_alpha[0] + 88, c_alpha[1])
    s.text(c_alpha[0] + 116, c_alpha[1], "R", 47)
    s.line(c_alpha[0], c_alpha[1] + 33, c_alpha[0], c_alpha[1] + 78)
    s.text(c_alpha[0], c_alpha[1] + 107, "H", 45)


def scheme_edman(s):
    y = 270
    draw_pitc(s, 180, y)
    plus(s, 620, y)
    draw_horizontal_aa(s, 900, y)
    arrow(s, 1190, 1290, y, label_top="弱碱中")
    draw_ptc(s, 1380, y)
    arrow(s, 2200, 2300, y, label_top="H⁺", label_bottom="(CH₃NO₂)")
    draw_pth(s, 2630, 210)
    caption(s, 280, 570, "异硫氰酸苯酯")
    caption(s, 1680, 570, "苯氨基硫甲酰氨基酸（PTC-氨基酸）", 38)
    caption(s, 2650, 570, "苯乙内酰硫脲氨基酸（PTH-氨基酸）", 38)


def draw_aldehyde(s, x, y):
    s.text(x, y, "C", 57)
    s.line(x, y - 36, x, y - 82)
    draw_r_prime(s, x, y - 113, 54)
    s.line(x, y + 36, x, y + 82)
    s.text(x, y + 113, "H", 52)
    double_bond(s, x + 36, y - 7, x + 84, y - 7, 10, 3.5)
    s.text(x + 118, y - 7, "O", 55)


def draw_schiff_product(s, x, y):
    s.text(x, y, "C", 56)
    s.line(x, y - 36, x, y - 82)
    draw_r_prime(s, x, y - 113, 53)
    s.line(x, y + 36, x, y + 82)
    s.text(x, y + 113, "H", 51)
    nx = x + 92
    double_bond(s, x + 36, y, nx - 32, y, 9, 3.5)
    s.text(nx, y, "N", 55)
    chx = nx + 108
    s.line(nx + 30, y, chx - 27, y)
    s.text(chx, y, "C", 56)
    s.line(chx, y - 31, chx, y - 74)
    draw_vertical_connected_label(s, chx, y - 104, "COOH", 52)
    s.line(chx, y + 31, chx, y + 74)
    s.text(chx, y + 104, "R", 52)
    s.line(chx + 27, y, chx + 76, y)
    s.text(chx + 82, y, "H", 50, anchor="start")


def scheme_schiff(s):
    y = 270
    draw_aldehyde(s, 240, y)
    plus(s, 450, y)
    draw_fischer_aa(s, 700, y)
    equilibrium_arrow(s, 895, 1015, y)
    draw_schiff_product(s, 1120, y)
    caption(s, 240, 470, "醛")
    caption(s, 700, 470, "氨基酸")
    caption(s, 1290, 470, "希夫碱")


SCHEMES = [
    ("01_氨基酸与亚硝酸反应", 1900, 520, scheme_nitrous),
    ("02_氨基酸酰基化反应", 1580, 550, scheme_acylation),
    ("03_DNFB与氨基酸反应", 2320, 550, scheme_dnfb),
    ("04_Edman降解反应", 3100, 650, scheme_edman),
    ("05_氨基酸形成希夫碱", 1480, 550, scheme_schiff),
]


def validate_structures():
    checks = {
        "丙氨酸": ("N[C@@H](C)C(=O)O", "C3H7NO2"),
        "亚硝酸": ("N(=O)O", "HNO2"),
        "乳酸": ("C[C@@H](O)C(=O)O", "C3H6O3"),
        "乙酰氯": ("CC(=O)Cl", "C2H3ClO"),
        "N-乙酰丙氨酸": ("CC(=O)N[C@@H](C)C(=O)O", "C5H9NO3"),
        "DNFB": ("Fc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1", "C6H3FN2O4"),
        "DNP-丙氨酸": ("O=[N+]([O-])c1ccc([N+](=O)[O-])c(N[C@@H](C)C(=O)O)c1", "C9H9N3O6"),
        "异硫氰酸苯酯": ("S=C=Nc1ccccc1", "C7H5NS"),
        "PTC-丙氨酸": ("O=C(O)[C@@H](C)NC(=S)Nc1ccccc1", "C10H12N2O2S"),
        "PTH-丙氨酸": ("O=C1[C@@H](C)NC(=S)N1c1ccccc1", "C10H10N2OS"),
        "乙醛": ("CC=O", "C2H4O"),
        "乙醛-丙氨酸希夫碱": ("CC=N[C@@H](C)C(=O)O", "C5H9NO2"),
    }
    results = []
    for name, (smiles, expected) in checks.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"RDKit无法解析：{name}")
        formula = rdMolDescriptors.CalcMolFormula(mol)
        if formula != expected:
            raise ValueError(f"{name}分子式校验失败：{formula} != {expected}")
        results.append((name, formula))
    return results


def export_png(svg_text, output_path, dpi):
    try:
        import cairosvg
        from PIL import Image
    except ImportError as exc:
        raise RuntimeError(
            "PNG export requires CairoSVG and Pillow. Install requirements.txt "
            "or run with --formats svg."
        ) from exc

    rendered = BytesIO()
    cairosvg.svg2png(bytestring=svg_text.encode("utf-8"), write_to=rendered)
    rendered.seek(0)
    with Image.open(rendered) as image:
        image.save(output_path, dpi=(dpi, dpi))


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Draw five textbook-ready amino-acid reaction schemes with "
            "validated structures and standardized typography."
        )
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUT,
        help="Directory for white and transparent SVG/PNG files.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("svg", "png"),
        default=("svg", "png"),
        help="Output formats. Defaults to both SVG and PNG.",
    )
    parser.add_argument(
        "--png-dpi",
        type=int,
        default=300,
        help="PNG resolution metadata. Defaults to 300 dpi.",
    )
    parser.add_argument(
        "--times-font",
        type=Path,
        default=CHEM_FONT_FILE,
        help="Times New Roman-compatible font file used for chemical text.",
    )
    parser.add_argument(
        "--simsun-font",
        type=Path,
        default=CN_FONT_FILE,
        help="SimSun-compatible font file used for Chinese captions.",
    )
    return parser.parse_args()


def main():
    global OUT, CHEM_FONT_FILE, CN_FONT_FILE
    args = parse_args()
    OUT = args.output_dir.resolve()
    OUT.mkdir(parents=True, exist_ok=True)
    CHEM_FONT_FILE = args.times_font.resolve()
    CN_FONT_FILE = args.simsun_font.resolve()
    for font_path in (CHEM_FONT_FILE, CN_FONT_FILE):
        if not font_path.is_file():
            raise FileNotFoundError(f"Font file not found: {font_path}")
    verified = validate_structures()
    for name, width, height, builder in SCHEMES:
        for white, suffix in ((True, "白底"), (False, "透明背景")):
            svg = SVG(width, height, white)
            builder(svg)
            svg_text = svg.finish()
            if "svg" in args.formats:
                svg_path = OUT / f"{name}_{suffix}.svg"
                svg_path.write_text(svg_text, encoding="utf-8")
                print(svg_path)
            if "png" in args.formats:
                png_path = OUT / f"{name}_{suffix}_{args.png_dpi}dpi.png"
                export_png(svg_text, png_path, args.png_dpi)
                print(png_path)
    print("RDKit结构校验：")
    for item in verified:
        print(item)


if __name__ == "__main__":
    main()
