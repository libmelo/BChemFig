---
name: amino-acid-fischer-projections
description: Draw or review textbook-ready Fischer projection diagrams of amino acids, including zwitterionic states, explicit heavy-atom connectivity, stereochemistry, Times New Roman chemical labels, Chinese captions, and grouped SVG/PNG output. Use for 氨基酸投影式、20种氨基酸结构组图 or corrections to such diagrams; do not use for ordinary skeletal formulas that are not Fischer projections.
---

# Amino Acid Fischer Projections

Produce chemically correct, visually consistent amino-acid Fischer projections suitable for Chinese biochemistry textbooks.

## Required procedure

1. Read [references/drawing_rules.md](references/drawing_rules.md) before drawing, modifying, or reviewing an amino-acid projection. Treat its stereochemical and connectivity rules as invariants.
2. Determine the requested amino acids, charge convention, labels, grouping, and output formats. If the user asks for the standard 20 without further detail, use the five groups and zwitterionic convention in the reference.
3. Prefer [scripts/draw_amino_acid_fischer.py](scripts/draw_amino_acid_fischer.py) for the canonical 20-amino-acid set. Adapt a copy for a different subset or layout; preserve the chemistry and typography invariants.
4. Export both SVG and 300-dpi PNG unless the user requests one format. Generate individual diagrams and the requested group composites. Put group composites in a clearly named, separate folder when individual files are also present.
5. Inspect every final group image at original resolution. Verify every heavy-atom bond, carbon label, charge, stereocenter, ring, double bond, caption, font, baseline, and crop. Do not report completion based only on successful execution.

## Deterministic rendering

The helper requires Python, RDKit, and Pillow. It searches `RDKIT_VENDOR_PATH` and ancestor `.vendor/rdkit` directories before using the environment installation.

```powershell
python scripts/draw_amino_acid_fischer.py --output-dir <output-directory>
```

On Windows it defaults to `C:\Windows\Fonts\times.ttf` and `C:\Windows\Fonts\simsun.ttc`. On other systems, pass compatible font files explicitly:

```powershell
python scripts/draw_amino_acid_fischer.py --output-dir <output-directory> --times-font <font-file> --simsun-font <font-file>
```

## Non-negotiable rendering invariants

- Use one text engine, one font file, one point size, and one baseline for the main `C` and `H` characters in every `CH`, `CH₂`, `CH₃`, or `H₃C` label. Draw numeric subscripts separately at a smaller size and lowered position. Merely assigning the same nominal size in two renderers is insufficient.
- Keep the carbon symbol adjacent to the bond. Use `H₃C—` when a terminal methyl bond leaves to the right; use `—CH₃` when it enters from the left or above as appropriate.
- Let RDKit establish atom-aware bond gaps. When replacing a carbon glyph for unified text rendering, mask only the original glyph area; a large white rectangle will create broken-looking bonds.
- Chemical symbols and Latin abbreviations use Times New Roman or the user-specified equivalent. Chinese captions use SimSun or the user-specified Chinese font. Do not fake subscripts with full-size baseline digits.
- Aromatic ring carbon vertices may remain implicit by standard skeletal convention. All other non-hydrogen atoms and all bonds between them must be visible. `COO⁻` may be abbreviated when the user permits it.

## Review behavior

When asked to diagnose an existing image, identify the exact chemical or typographic cause before editing. Typical causes include a bond ending on `H` instead of `C`, a carbon hidden by an abbreviation, mixed rendering engines, imidazole attachment at the wrong carbon, or an overlarge glyph mask. Correct the source and regenerate every affected individual and group file.
