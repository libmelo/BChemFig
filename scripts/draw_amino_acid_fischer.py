from pathlib import Path
import argparse
import io
import math
import os
import sys


def find_rdkit_vendor():
    """Find a locally vendored RDKit without tying the skill to one workspace."""
    candidates = []
    configured = os.environ.get("RDKIT_VENDOR_PATH")
    if configured:
        candidates.append(Path(configured))
    for start in (Path.cwd(), Path(__file__).resolve().parent):
        candidates.extend(parent / ".vendor" / "rdkit" for parent in (start, *start.parents))
    for candidate in candidates:
        if candidate.is_dir():
            return candidate
    return None


VENDOR = find_rdkit_vendor()
if VENDOR is not None:
    sys.path.insert(0, str(VENDOR))
    if hasattr(os, "add_dll_directory"):
        for dll_dir in (VENDOR / "bin", VENDOR / "rdkit.libs"):
            if dll_dir.is_dir():
                os.add_dll_directory(str(dll_dir))

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageChops, ImageDraw, ImageFont


OUT = Path.cwd() / "amino_acid_fischer_output"
TIMES = r"C:\Windows\Fonts\times.ttf"
SIMSUN = r"C:\Windows\Fonts\simsun.ttc"
ATOM_TEXT_SIZE = 44
ATOM_SUBSCRIPT_SIZE = 27


class Diagram:
    def __init__(self):
        self.rw = Chem.RWMol()
        self.xy = {}

    def add(self, label, x, y, atomic_num=0):
        atom = Chem.Atom(atomic_num)
        if label is not None:
            atom.SetProp("atomLabel", label)
        idx = self.rw.AddAtom(atom)
        self.xy[idx] = (x, y)
        return idx

    def bond(self, a, b, order=Chem.BondType.SINGLE):
        self.rw.AddBond(a, b, order)

    def add_element(
        self, symbol, x, y, hydrogens=0, charge=0, force_label=False
    ):
        atom = Chem.Atom(symbol)
        atom.SetFormalCharge(charge)
        atom.SetNumExplicitHs(hydrogens)
        atom.SetNoImplicit(True)
        if force_label:
            atom.SetProp("atomLabel", symbol)
        if symbol == "C" and hydrogens:
            atom.SetIntProp("_TextHydrogenCount", hydrogens)
        idx = self.rw.AddAtom(atom)
        self.xy[idx] = (x, y)
        return idx

    def finish(self):
        mol = self.rw.GetMol()
        mol.UpdatePropertyCache(strict=False)
        conf = Chem.Conformer(mol.GetNumAtoms())
        for idx, (x, y) in self.xy.items():
            conf.SetAtomPosition(idx, (x, y, 0.0))
        mol.RemoveAllConformers()
        mol.AddConformer(conf)
        return mol


def core():
    d = Diagram()
    ca = d.add_element("C", 0.0, 0.0, force_label=True)
    carboxylate = d.add("COO<sup>-</sup>", 0.0, 0.82)
    amino = d.add_element("N", -0.98, 0.0, hydrogens=3, charge=1)
    hydrogen = d.add("H", 0.92, 0.0)
    d.bond(ca, carboxylate)
    d.bond(ca, amino)
    d.bond(ca, hydrogen)
    return d, ca


def node(d, parent, label, x, y, order=Chem.BondType.SINGLE):
    carbon_hydrogens = {
        "C": 0,
        "CH": 1,
        "CH<sub>2</sub>": 2,
        "CH<sub>3</sub>": 3,
    }
    if label in carbon_hydrogens:
        idx = d.add_element(
            "C", x, y, hydrogens=carbon_hydrogens[label], force_label=True
        )
    else:
        idx = d.add(label, x, y)
    d.bond(parent, idx, order)
    return idx


def element_node(
    d, parent, symbol, x, y, hydrogens=0, charge=0,
    order=Chem.BondType.SINGLE,
):
    idx = d.add_element(symbol, x, y, hydrogens=hydrogens, charge=charge)
    d.bond(parent, idx, order)
    return idx


def vertical_chain(d, parent, labels, start_y=-0.82, step=0.78):
    current = parent
    y = start_y
    nodes = []
    for label in labels:
        current = node(d, current, label, 0.0, y)
        nodes.append(current)
        y -= step
    return nodes, y


def terminal_carboxylate(d, parent, carbon_y):
    carboxylate = node(d, parent, "COO<sup>-</sup>", 0.0, carbon_y)
    return carboxylate


def terminal_amide(d, parent, carbon_y):
    c = node(d, parent, "C", 0.0, carbon_y)
    o = element_node(d, c, "O", -0.65, carbon_y - 0.62, order=Chem.BondType.DOUBLE)
    n = element_node(d, c, "N", 0.72, carbon_y - 0.62, hydrogens=2)
    return c, o, n


def append_ring_fragment(
    d, parent, smiles, attach_atom, target, scale=0.48, mirror_x=False
):
    frag = Chem.MolFromSmiles(smiles)
    if frag is None:
        raise RuntimeError(f"Cannot build fragment: {smiles}")
    rdDepictor.Compute2DCoords(frag)
    Chem.Kekulize(frag, clearAromaticFlags=True)
    conf = frag.GetConformer()
    coords = [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in range(frag.GetNumAtoms())]
    ax, ay = coords[attach_atom]
    cx = sum(x for x, _ in coords) / len(coords)
    cy = sum(y for _, y in coords) / len(coords)
    current_angle = math.atan2(cy - ay, cx - ax)
    rotation = -math.pi / 2 - current_angle
    cos_r, sin_r = math.cos(rotation), math.sin(rotation)

    mapping = {}
    for atom in frag.GetAtoms():
        old = atom.GetIdx()
        x, y = coords[old]
        dx, dy = x - ax, y - ay
        rx = dx * cos_r - dy * sin_r
        ry = dx * sin_r + dy * cos_r
        if mirror_x:
            rx = -rx
        tx = target[0] + scale * rx
        ty = target[1] + scale * ry
        new_atom = Chem.Atom(atom)
        mapping[old] = d.rw.AddAtom(new_atom)
        d.xy[mapping[old]] = (tx, ty)
    for bond in frag.GetBonds():
        d.rw.AddBond(
            mapping[bond.GetBeginAtomIdx()],
            mapping[bond.GetEndAtomIdx()],
            bond.GetBondType(),
        )
    d.bond(parent, mapping[attach_atom])
    return mapping


def simple_vertical(labels):
    d, ca = core()
    vertical_chain(d, ca, labels)
    return d.finish()


def glycine():
    return simple_vertical(["H"])


def alanine():
    return simple_vertical(["CH<sub>3</sub>"])


def valine():
    d, ca = core()
    beta = node(d, ca, "CH", 0.0, -0.82)
    node(d, beta, "CH<sub>3</sub>", -0.72, -1.48)
    node(d, beta, "CH<sub>3</sub>", 0.72, -1.48)
    return d.finish()


def leucine():
    d, ca = core()
    cb = node(d, ca, "CH<sub>2</sub>", 0.0, -0.82)
    cg = node(d, cb, "CH", 0.0, -1.60)
    node(d, cg, "CH<sub>3</sub>", -0.72, -2.26)
    node(d, cg, "CH<sub>3</sub>", 0.72, -2.26)
    return d.finish()


def isoleucine():
    d, ca = core()
    cb = node(d, ca, "C", 0.0, -0.82)
    node(d, cb, "CH<sub>3</sub>", -0.88, -0.82)
    node(d, cb, "H", 0.86, -0.82)
    cg = node(d, cb, "CH<sub>2</sub>", 0.0, -1.64)
    node(d, cg, "CH<sub>3</sub>", 0.0, -2.42)
    return d.finish()


def aspartate():
    d, ca = core()
    cb = node(d, ca, "CH<sub>2</sub>", 0.0, -0.82)
    terminal_carboxylate(d, cb, -1.60)
    return d.finish()


def glutamate():
    d, ca = core()
    nodes, y = vertical_chain(d, ca, ["CH<sub>2</sub>", "CH<sub>2</sub>"])
    terminal_carboxylate(d, nodes[-1], y)
    return d.finish()


def lysine():
    d, ca = core()
    nodes, y = vertical_chain(
        d, ca,
        ["CH<sub>2</sub>", "CH<sub>2</sub>", "CH<sub>2</sub>", "CH<sub>2</sub>"],
    )
    element_node(d, nodes[-1], "N", 0.0, y, hydrogens=3, charge=1)
    return d.finish()


def arginine():
    d, ca = core()
    nodes, y = vertical_chain(d, ca, ["CH<sub>2</sub>", "CH<sub>2</sub>", "CH<sub>2</sub>"])
    chain_n = element_node(d, nodes[-1], "N", 0.0, y, hydrogens=1)
    guan_c = node(d, chain_n, "C", 0.0, y - 0.78)
    element_node(
        d, guan_c, "N", -0.78, y - 1.42,
        hydrogens=2, charge=1, order=Chem.BondType.DOUBLE,
    )
    element_node(d, guan_c, "N", 0.78, y - 1.42, hydrogens=2)
    return d.finish()


def serine():
    d, ca = core()
    cb = node(d, ca, "CH<sub>2</sub>", 0.0, -0.82)
    element_node(d, cb, "O", 0.0, -1.60, hydrogens=1)
    return d.finish()


def threonine():
    d, ca = core()
    cb = node(d, ca, "C", 0.0, -0.82)
    node(d, cb, "H", -0.86, -0.82)
    element_node(d, cb, "O", 0.86, -0.82, hydrogens=1)
    node(d, cb, "CH<sub>3</sub>", 0.0, -1.64)
    return d.finish()


def cysteine():
    d, ca = core()
    cb = node(d, ca, "CH<sub>2</sub>", 0.0, -0.82)
    element_node(d, cb, "S", 0.0, -1.60, hydrogens=1)
    return d.finish()


def methionine():
    d, ca = core()
    nodes, y = vertical_chain(d, ca, ["CH<sub>2</sub>", "CH<sub>2</sub>"])
    sulfur = element_node(d, nodes[-1], "S", 0.0, y)
    node(d, sulfur, "CH<sub>3</sub>", 0.0, y - 0.78)
    return d.finish()


def asparagine():
    d, ca = core()
    cb = node(d, ca, "CH<sub>2</sub>", 0.0, -0.82)
    terminal_amide(d, cb, -1.60)
    return d.finish()


def glutamine():
    d, ca = core()
    nodes, y = vertical_chain(d, ca, ["CH<sub>2</sub>", "CH<sub>2</sub>"])
    terminal_amide(d, nodes[-1], y)
    return d.finish()


def phenylalanine():
    d, ca = core()
    cb = node(d, ca, "CH<sub>2</sub>", 0.0, -0.82)
    append_ring_fragment(d, cb, "c1ccccc1", 0, (0.0, -1.57), 0.48)
    return d.finish()


def tryptophan():
    d, ca = core()
    cb = node(d, ca, "CH<sub>2</sub>", 0.0, -0.82)
    append_ring_fragment(
        d, cb, "c1c[nH]c2ccccc12", 0, (0.0, -1.57),
        0.46, mirror_x=True,
    )
    return d.finish()


def tyrosine():
    d, ca = core()
    cb = node(d, ca, "CH<sub>2</sub>", 0.0, -0.82)
    mapping = append_ring_fragment(d, cb, "c1ccccc1", 0, (0.0, -1.57), 0.48)
    para = mapping[3]
    px, py = d.xy[para]
    cx = sum(d.xy[mapping[i]][0] for i in mapping) / len(mapping)
    cy = sum(d.xy[mapping[i]][1] for i in mapping) / len(mapping)
    vx, vy = px - cx, py - cy
    length = math.hypot(vx, vy)
    oh = d.add_element(
        "O", px + 0.62 * vx / length, py + 0.62 * vy / length,
        hydrogens=1,
    )
    d.bond(para, oh)
    return d.finish()


def histidine():
    d, ca = core()
    cb = node(d, ca, "CH<sub>2</sub>", 0.0, -0.82)
    append_ring_fragment(d, cb, "c1cnc[nH]1", 0, (0.0, -1.57), 0.50)
    return d.finish()


def proline():
    d = Diagram()
    ca = d.add_element("C", 0.0, 0.0, force_label=True)
    carboxylate = d.add("COO<sup>-</sup>", 0.0, 0.82)
    h = d.add("H", 0.92, 0.0)
    n = d.add_element("N", -0.98, 0.0, hydrogens=2, charge=1)
    cb = d.add_element("C", 0.0, -0.82, hydrogens=2, force_label=True)
    cg = d.add_element("C", -0.68, -1.42, hydrogens=2, force_label=True)
    cd = d.add_element("C", -1.26, -0.70, hydrogens=2, force_label=True)
    d.bond(ca, carboxylate)
    d.bond(ca, h)
    d.bond(ca, n)
    d.bond(ca, cb)
    d.bond(cb, cg)
    d.bond(cg, cd)
    d.bond(cd, n)
    return d.finish()


def caption(name, three, one):
    return [(name, "simsun"), (f" ({three}, {one})", "times")]


ALL_GROUPS = [
    ("01_Gly_Ala_Val_Leu_Ile", [
        ("glycine_zwitterionic_fischer", glycine(), caption("甘氨酸", "Gly", "G")),
        ("alanine_zwitterionic_fischer", alanine(), caption("丙氨酸", "Ala", "A")),
        ("valine_zwitterionic_fischer", valine(), caption("缬氨酸", "Val", "V")),
        ("leucine_zwitterionic_fischer", leucine(), caption("亮氨酸", "Leu", "L")),
        ("isoleucine_zwitterionic_fischer", isoleucine(), caption("异亮氨酸", "Ile", "I")),
    ]),
    ("02_Asp_Glu_Lys_Arg_Ser_Thr", [
        ("aspartate_zwitterionic_fischer", aspartate(), caption("天冬氨酸", "Asp", "D")),
        ("glutamate_zwitterionic_fischer", glutamate(), caption("谷氨酸", "Glu", "E")),
        ("lysine_zwitterionic_fischer", lysine(), caption("赖氨酸", "Lys", "K")),
        ("arginine_zwitterionic_fischer", arginine(), caption("精氨酸", "Arg", "R")),
        ("serine_zwitterionic_fischer", serine(), caption("丝氨酸", "Ser", "S")),
        ("threonine_zwitterionic_fischer", threonine(), caption("苏氨酸", "Thr", "T")),
    ]),
    ("03_Cys_Met_Asn_Gln", [
        ("cysteine_zwitterionic_fischer", cysteine(), caption("半胱氨酸", "Cys", "C")),
        ("methionine_zwitterionic_fischer", methionine(), caption("甲硫氨酸", "Met", "M")),
        ("asparagine_zwitterionic_fischer", asparagine(), caption("天冬酰胺", "Asn", "N")),
        ("glutamine_zwitterionic_fischer", glutamine(), caption("谷氨酰胺", "Gln", "Q")),
    ]),
    ("04_Phe_Trp_Tyr", [
        ("phenylalanine_zwitterionic_fischer", phenylalanine(), caption("苯丙氨酸", "Phe", "F")),
        ("tryptophan_zwitterionic_fischer", tryptophan(), caption("色氨酸", "Trp", "W")),
        ("tyrosine_zwitterionic_fischer", tyrosine(), caption("酪氨酸", "Tyr", "Y")),
    ]),
    ("05_His_Pro", [
        ("histidine_zwitterionic_fischer", histidine(), caption("组氨酸", "His", "H")),
        ("proline_zwitterionic_fischer", proline(), caption("脯氨酸", "Pro", "P")),
    ]),
]


def configure(options):
    options.bondLineWidth = 2
    options.padding = 0.20
    options.baseFontSize = 0.26
    options.minFontSize = 38
    options.maxFontSize = 50
    options.fixedBondLength = 92
    options.addStereoAnnotation = False
    options.explicitMethyl = True
    options.useBWAtomPalette()
    options.fontFile = TIMES


def collect_carbon_hydrogen_overlays(mol, drawer):
    overlays = []
    for atom in mol.GetAtoms():
        if not atom.HasProp("_TextHydrogenCount"):
            continue
        count = atom.GetIntProp("_TextHydrogenCount")
        position = drawer.GetDrawCoords(atom.GetIdx())
        neighbours = list(atom.GetNeighbors())
        prefix = False
        # For a terminal carbon whose heavy-atom neighbour lies to its right,
        # write H3C so that the C remains adjacent to the bond. Otherwise use
        # CHn. Internal carbons normally use the right-hand suffix form.
        if len(neighbours) == 1:
            neighbour_position = drawer.GetDrawCoords(neighbours[0].GetIdx())
            prefix = neighbour_position.x > position.x + 2
        overlays.append((position.x, position.y, count, prefix))
    return overlays


def label_geometry(x, y, count, prefix, main_font, subscript_font):
    """Place C and H on one shared baseline using one set of font metrics."""
    c_width = main_font.getlength("C")
    h_width = main_font.getlength("H")
    digit_width = subscript_font.getlength(str(count)) if count > 1 else 0
    c_box = main_font.getbbox("C", anchor="ms")
    baseline_y = y - (c_box[1] + c_box[3]) / 2

    if prefix:
        digit_x = x - c_width / 2 - digit_width / 2
        h_x = x - c_width / 2 - digit_width - h_width / 2
    else:
        h_x = x + c_width / 2 + h_width / 2
        digit_x = x + c_width / 2 + h_width + digit_width / 2
    return x, h_x, digit_x, baseline_y


def add_svg_carbon_hydrogen_overlays(svg_text, overlays):
    main_font = ImageFont.truetype(TIMES, ATOM_TEXT_SIZE)
    subscript_font = ImageFont.truetype(TIMES, ATOM_SUBSCRIPT_SIZE)
    parts = []
    for x, y, count, prefix in overlays:
        c_x, h_x, digit_x, baseline_y = label_geometry(
            x, y, count, prefix, main_font, subscript_font
        )
        # Hide RDKit's C glyph, but retain the bond layout that it established.
        # C and H are then redrawn together by the same text engine.
        parts.append(
            f'<rect x="{x - 23:.1f}" y="{y - 22:.1f}" width="46" height="44" '
            f'fill="white"/>'
        )
        parts.append(
            f'<text x="{c_x:.1f}" y="{baseline_y:.1f}" text-anchor="middle" '
            f'font-family="Times New Roman" font-size="{ATOM_TEXT_SIZE}">C</text>'
        )
        parts.append(
            f'<text x="{h_x:.1f}" y="{baseline_y:.1f}" text-anchor="middle" '
            f'font-family="Times New Roman" font-size="{ATOM_TEXT_SIZE}">H</text>'
        )
        if count > 1:
            parts.append(
                f'<text x="{digit_x:.1f}" y="{baseline_y + 11:.1f}" '
                f'text-anchor="middle" '
                f'font-family="Times New Roman" font-size="{ATOM_SUBSCRIPT_SIZE}">{count}</text>'
            )
    return svg_text.replace("</svg>", "\n".join(parts) + "</svg>")


def add_png_carbon_hydrogen_overlays(png_bytes, overlays):
    image = Image.open(io.BytesIO(png_bytes)).convert("RGB")
    draw = ImageDraw.Draw(image)
    main_font = ImageFont.truetype(TIMES, ATOM_TEXT_SIZE)
    subscript_font = ImageFont.truetype(TIMES, ATOM_SUBSCRIPT_SIZE)
    for x, y, count, prefix in overlays:
        c_x, h_x, digit_x, baseline_y = label_geometry(
            x, y, count, prefix, main_font, subscript_font
        )
        draw.rectangle((x - 23, y - 22, x + 23, y + 22), fill="white")
        draw.text((c_x, baseline_y), "C", font=main_font, fill="black", anchor="ms")
        draw.text((h_x, baseline_y), "H", font=main_font, fill="black", anchor="ms")
        if count > 1:
            draw.text(
                (digit_x, baseline_y + 11), str(count), font=subscript_font,
                fill="black", anchor="ms",
            )
    return image


def nonwhite_bbox(image, margin=20):
    rgb = image.convert("RGB")
    diff = ImageChops.difference(rgb, Image.new("RGB", rgb.size, "white"))
    box = diff.getbbox() or (0, 0, rgb.width, rgb.height)
    l, t, r, b = box
    return max(0, l-margin), max(0, t-margin), min(rgb.width, r+margin), min(rgb.height, b+margin)


def svg_inner(text):
    start = text.find(">", text.find("<svg")) + 1
    return text[start:text.rfind("</svg>")]


def draw_mixed(draw, center_x, baseline_y, segments, size):
    fonts = {
        "times": ImageFont.truetype(TIMES, size),
        "simsun": ImageFont.truetype(SIMSUN, size),
    }
    widths = [draw.textlength(text, font=fonts[kind]) for text, kind in segments]
    x = center_x - sum(widths)/2
    for (text, kind), width in zip(segments, widths):
        draw.text((x, baseline_y), text, font=fonts[kind], fill="black", anchor="ls")
        x += width


def render_all():
    for group_id, specs in ALL_GROUPS:
        exported = []
        for stem, mol, segments in specs:
            w, h = 1400, 1100
            sd = rdMolDraw2D.MolDraw2DSVG(w, h)
            configure(sd.drawOptions())
            sd.DrawMolecule(mol)
            carbon_hydrogen_overlays = collect_carbon_hydrogen_overlays(mol, sd)
            sd.FinishDrawing()
            svg_path = OUT / f"{stem}.svg"
            svg_path.write_text(
                add_svg_carbon_hydrogen_overlays(
                    sd.GetDrawingText(), carbon_hydrogen_overlays
                ),
                encoding="utf-8",
            )

            pd = rdMolDraw2D.MolDraw2DCairo(w, h)
            configure(pd.drawOptions())
            pd.DrawMolecule(mol)
            carbon_hydrogen_overlays = collect_carbon_hydrogen_overlays(mol, pd)
            pd.FinishDrawing()
            png_path = OUT / f"{stem}.png"
            add_png_carbon_hydrogen_overlays(
                pd.GetDrawingText(), carbon_hydrogen_overlays
            ).save(png_path, dpi=(300, 300))
            exported.append((svg_path, png_path, segments))

        count = len(exported)
        spacing = 760
        group_w, group_h = count * spacing, 980
        centers = [spacing/2 + i*spacing for i in range(count)]
        sources = [Image.open(p).convert("RGB") for _, p, _ in exported]
        boxes = [nonwhite_bbox(im, 22) for im in sources]
        max_w = max(b[2]-b[0] for b in boxes)
        max_h = max(b[3]-b[1] for b in boxes)
        scale = min(1.08, (spacing-70)/max_w, 690/max_h)
        structure_y = 405
        canvas = Image.new("RGB", (group_w, group_h), "white")
        draw = ImageDraw.Draw(canvas)
        for cx, source, box, (_, _, segments) in zip(centers, sources, boxes, exported):
            crop = source.crop(box)
            panel = crop.resize((round(crop.width*scale), round(crop.height*scale)), Image.Resampling.LANCZOS)
            canvas.paste(panel, (round(cx-panel.width/2), round(structure_y-panel.height/2)))
            draw_mixed(draw, cx, 900, segments, 35)
        group_png = OUT / f"amino_acids_group_{group_id}.png"
        canvas.save(group_png, dpi=(300, 300))

        panels = []
        for cx, box, (svg_path, _, _) in zip(centers, boxes, exported):
            sw, sh = box[2]-box[0], box[3]-box[1]
            pw, ph = round(sw*scale), round(sh*scale)
            panels.append(
                f'<svg x="{cx-pw/2}" y="{structure_y-ph/2}" width="{pw}" height="{ph}" '
                f'viewBox="{box[0]} {box[1]} {sw} {sh}">{svg_inner(svg_path.read_text(encoding="utf-8"))}</svg>'
            )
        captions = "\n".join(
            f'<text x="{cx}" y="900" text-anchor="middle" font-family="SimSun, 宋体" font-size="35">'
            f'{segments[0][0]}<tspan font-family="Times New Roman">{segments[1][0]}</tspan></text>'
            for cx, (_, _, segments) in zip(centers, exported)
        )
        group_svg = OUT / f"amino_acids_group_{group_id}.svg"
        group_svg.write_text(
            f'<?xml version="1.0" encoding="UTF-8"?><svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{group_w}px" height="{group_h}px" viewBox="0 0 {group_w} {group_h}">'
            f'<rect width="100%" height="100%" fill="white"/>{"".join(panels)}{captions}</svg>',
            encoding="utf-8",
        )
        print(group_svg)
        print(group_png)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Draw textbook-ready zwitterionic Fischer projections for the 20 standard amino acids."
    )
    parser.add_argument(
        "--output-dir", type=Path, default=OUT,
        help="Directory for individual and five grouped PNG/SVG files.",
    )
    parser.add_argument(
        "--times-font", default=TIMES,
        help="Times New Roman-compatible font file used for chemical text.",
    )
    parser.add_argument(
        "--simsun-font", default=SIMSUN,
        help="SimSun-compatible font file used for Chinese captions.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    OUT = args.output_dir.resolve()
    OUT.mkdir(parents=True, exist_ok=True)
    TIMES = str(Path(args.times_font).resolve())
    SIMSUN = str(Path(args.simsun_font).resolve())
    for font_path in (TIMES, SIMSUN):
        if not Path(font_path).is_file():
            raise FileNotFoundError(f"Font file not found: {font_path}")
    render_all()
