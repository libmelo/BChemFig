---
name: bchemfig
description: Create or review textbook-ready biochemical figures for Chinese biochemistry and food-biochemistry textbooks, including chemical structures and reactions, data curves, separation apparatus and principles, protein structural models, sequence/topology comparisons, workflows, tables, and Word-editable equations. Use for 生化教材插图、蛋白质或氨基酸图、实验原理图及其科学性和版式修订; do not use for decorative biological artwork or unrelated general-purpose graphics.
metadata:
  short-description: 制作与审校生化教材插图
---

# BChemFig

Produce scientifically correct, visually restrained, print-ready biochemical teaching figures. Treat a supplied textbook image as a content reference, not as an automatically correct structure or layout.

## Route to the relevant rules

Read [references/common_standards.md](references/common_standards.md) for every task, then read only the references needed for the requested artifact:

- Chemical structures, Fischer projections, amino-acid plates: [references/drawing_rules.md](references/drawing_rules.md).
- Chemical reaction schemes: [references/reaction_scheme_rules.md](references/reaction_scheme_rules.md).
- Curves, quantitative plots, tables, or Word equations: [references/data_plots_tables_equations.md](references/data_plots_tables_equations.md).
- Separation principles, experimental apparatus, or teaching flowcharts: [references/apparatus_and_flowcharts.md](references/apparatus_and_flowcharts.md).
- Protein secondary, supersecondary, tertiary, or molecular-surface models: [references/protein_structure_models.md](references/protein_structure_models.md).
- Residue sequences, mutation comparisons, disulfide topology, or short-peptide diagrams: [references/sequence_and_topology.md](references/sequence_and_topology.md).

If one plate combines several modes, apply every relevant reference. Scientific correctness overrides visual similarity to the supplied reference.

## Required workflow

1. Identify the teaching point, biochemical entities, protonation or experimental conditions, required labels, page context, and output formats. Resolve substantive ambiguity from authoritative data or ask only when different choices would materially change the figure.
2. Verify the scientific content before layout. Check molecular connectivity, stereochemistry, charge, reaction conditions, sequence numbering, structural parameters, data values, units, and process state as applicable.
3. Choose a professional method matched to the artifact:
   - RDKit or another chemical toolkit for concrete molecular graphs; manual atom-aware vector layout for generalized `R`/`R′` schemes.
   - Reproducible numerical code and SVG/plotting tools for curves and calibrated graphs.
   - PDB experimental structures and PyMOL for three-dimensional proteins.
   - Vector drawing for apparatus, mechanisms, topology, and workflows.
   - Native Word equations and tables when editability in Word is required.
4. Build a restrained black, white, and gray teaching graphic. Preserve meaningful hierarchy and compact spacing; do not decorate, over-label, or reproduce avoidable defects from the reference.
5. Export SVG plus 300-dpi PNG by default. Keep a white-background version and, when useful for Word placement, an equivalent transparent-background version. Provide a Word file for editable equations or tables when requested.
6. Inspect every final artifact at original resolution. Verify science, fonts, line attachment, arrow direction, labels, overlap, crop, and consistency across white and transparent variants. Fix the source and regenerate; do not patch only the final PNG.

## Reusable deterministic scripts

For the canonical 20-amino-acid Fischer set:

```powershell
python scripts/draw_amino_acid_fischer.py --output-dir <output-directory>
```

For the validated five-scheme amino-acid reaction set:

```powershell
python scripts/draw_amino_acid_reactions.py --output-dir <output-directory>
```

Both helpers search `RDKIT_VENDOR_PATH` and ancestor `.vendor/rdkit` directories. On Windows they default to Times New Roman and SimSun system fonts. Adapt a copy for a new subset or layout rather than weakening the chemical, typography, or attachment invariants in the maintained scripts.

## Completion standard

Do not report completion merely because a file rendered. A final figure must be chemically or biochemically correct, legible at intended textbook size, free of unintended overlap, visually consistent with the series, and delivered in the promised formats.
