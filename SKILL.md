---
name: amino-acid-fischer-projections
description: Draw or review textbook-ready amino-acid Fischer projections and amino-acid reaction schemes, with explicit heavy-atom connectivity, stereochemical checks, standardized arrows, Times New Roman chemical labels, Chinese captions, and SVG/PNG output. Use for 氨基酸投影式、20种氨基酸结构组图、氨基酸反应式 or corrections to these textbook figures; do not use for unrelated general-purpose chemical illustration.
---

# Amino Acid Textbook Figures

Produce chemically correct, visually consistent amino-acid structural figures suitable for Chinese biochemistry textbooks.

## Choose the drawing mode

- For a Fischer projection, a standard-20 set, or a grouped amino-acid structure plate, read [references/drawing_rules.md](references/drawing_rules.md) completely before acting.
- For an amino-acid reaction scheme, reaction-arrow correction, or textbook reaction plate, read [references/reaction_scheme_rules.md](references/reaction_scheme_rules.md) completely before acting.
- If one artifact mixes projections and reaction schemes, apply both references. Shared rules do not override reaction-specific chemistry or Fischer stereochemistry.

## Required procedure

1. Determine the requested structures, protonation convention, reaction conditions, labels, grouping, and output formats. Preserve the user's reference convention unless it is chemically wrong; explain and correct substantive errors.
2. Prefer the bundled deterministic script for a canonical figure. Adapt a copy for a different subset or layout while preserving the relevant chemistry, attachment, typography, and spacing invariants.
3. Use RDKit or another professional chemical toolkit to validate structures, formulae, aromaticity, valence, and stereochemistry where a concrete molecular graph is available. Generalized `R`/`R′` schemes still require manual atom-by-atom review.
4. Export SVG and 300-dpi PNG unless the user requests one format. Preserve a white-background version; add a transparent-background version when useful for Word or page layout.
5. Inspect every final image at original resolution. Verify every heavy-atom bond, attachment atom, carbon label, charge, stereocenter, ring, double bond, reaction arrow, condition label, caption, font, baseline, and crop. Do not report completion based only on successful execution.

## Deterministic rendering: Fischer projections

The helper requires Python, RDKit, and Pillow. It searches `RDKIT_VENDOR_PATH` and ancestor `.vendor/rdkit` directories before using the environment installation.

```powershell
python scripts/draw_amino_acid_fischer.py --output-dir <output-directory>
```

On Windows it defaults to `C:\Windows\Fonts\times.ttf` and `C:\Windows\Fonts\simsun.ttc`. On other systems, pass compatible font files explicitly:

```powershell
python scripts/draw_amino_acid_fischer.py --output-dir <output-directory> --times-font <font-file> --simsun-font <font-file>
```

## Deterministic rendering: reaction schemes

The reaction helper generates the validated five-scheme teaching set: nitrous-acid deamination, N-acylation, DNFB labeling, Edman degradation, and Schiff-base formation.

```powershell
python scripts/draw_amino_acid_reactions.py --output-dir <output-directory>
```

It searches `RDKIT_VENDOR_PATH` and ancestor `.vendor/rdkit` directories before using the environment installation. The default output includes SVG and PNG; PNG rendering requires CairoSVG from `requirements.txt`. Use `--formats svg` when a vector-only intermediate is desired. On non-Windows systems, pass Times New Roman-compatible and SimSun-compatible files with `--times-font` and `--simsun-font`.

## Non-negotiable rendering invariants

- Use one text engine, one font file, one point size, and one baseline for the main `C` and `H` characters in every `CH`, `CH₂`, `CH₃`, or `H₃C` label. Draw numeric subscripts separately at a smaller size and lowered position. Merely assigning the same nominal size in two renderers is insufficient.
- Keep the carbon symbol adjacent to the bond. Use `H₃C—` when a terminal methyl bond leaves to the right; use `—CH₃` when it enters from the left or above as appropriate.
- Let RDKit establish atom-aware bond gaps. When replacing a carbon glyph for unified text rendering, mask only the original glyph area; a large white rectangle will create broken-looking bonds.
- Chemical symbols and Latin abbreviations use Times New Roman or the user-specified equivalent. Chinese captions use SimSun or the user-specified Chinese font. Do not fake subscripts with full-size baseline digits.
- Aromatic ring carbon vertices may remain implicit by standard skeletal convention. All other non-hydrogen atoms and all bonds between them must be visible. `COO⁻` may be abbreviated when the user permits it.

## Review behavior

When asked to diagnose an existing image, identify the exact chemical or typographic cause before editing. Typical causes include a bond ending on `H` instead of `C`, a carbon hidden by an abbreviation, mixed rendering engines, a mispositioned `R′`, imidazole attachment at the wrong carbon, a distorted PTH ring, or an overlarge glyph mask. Correct the source and regenerate every affected individual, group, white-background, and transparent-background file.
