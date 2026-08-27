# BChemFig：氨基酸教材插图

一个面向中文生物化学与食品生物化学教材的 Codex skill，用于绘制和审查氨基酸 Fischer 投影式及典型反应式。项目同时提供化学规则、排版规则、可重复生成脚本和逐图质检清单。

## 主要能力

- 遵循 L-α-氨基酸 Fischer 投影的立体化学约定。
- 显式绘制非环重原子及其化学键，并核对双键、芳香环与杂环连接位点。
- 统一 `C`、`H` 的 Times New Roman 字体、字面尺寸和基线。
- 中文名称使用宋体，英文缩写使用 Times New Roman。
- 内置20种氨基酸的默认五组组图版式。
- 内置亚硝酸反应、N-酰基化、DNFB 标记、Edman 降解和希夫碱形成五类反应式。
- 反应式采用紧凑短箭头，严格保证键连接到正确的 C、N、O、S 原子。
- 对 `R′` 使用斜体 R 与独立上标撇号，避免不同公式中的基线漂移。
- 同时提供化学规则、排版规则、常见错误和逐图质检清单。

## 安装

将仓库克隆到 Codex 的个人 skills 目录：

```powershell
git clone git@github.com:libmelo/BChemFig.git "$env:USERPROFILE\.codex\skills\amino-acid-fischer-projections"
```

技能名称为：

```text
$amino-acid-fischer-projections
```

## Python 依赖

```powershell
python -m pip install -r requirements.txt
```

Windows 默认使用：

- `C:\Windows\Fonts\times.ttf`
- `C:\Windows\Fonts\simsun.ttc`

其他系统可通过命令行参数指定兼容字体文件。

## 直接生成20种氨基酸图

```powershell
python scripts/draw_amino_acid_fischer.py --output-dir output
```

指定字体：

```powershell
python scripts/draw_amino_acid_fischer.py \
  --output-dir output \
  --times-font <times-compatible-font-file> \
  --simsun-font <chinese-font-file>
```

脚本会生成20个单图及五组组图，每项均含 SVG 和 PNG。

## 直接生成氨基酸典型反应式

```powershell
python scripts/draw_amino_acid_reactions.py --output-dir output-reactions
```

默认生成五种反应式的白底、透明背景 SVG 和 300 dpi PNG。仅生成矢量文件时：

```powershell
python scripts/draw_amino_acid_reactions.py --output-dir output-reactions --formats svg
```

在非 Windows 系统上，可用 `--times-font` 与 `--simsun-font` 指定兼容字体文件。

## 技能结构

```text
SKILL.md
agents/openai.yaml
references/drawing_rules.md
references/reaction_scheme_rules.md
scripts/draw_amino_acid_fischer.py
scripts/draw_amino_acid_reactions.py
```

投影式规则见 [`references/drawing_rules.md`](references/drawing_rules.md)，反应式规则见 [`references/reaction_scheme_rules.md`](references/reaction_scheme_rules.md)。
