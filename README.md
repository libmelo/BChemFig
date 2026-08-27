# Amino Acid Fischer Projections

一个面向中文生物化学教材的 Codex skill，用于绘制和审查氨基酸 Fischer 投影式。默认生成20种标准 α-氨基酸的两性离子形式，同时输出单图和五组组图的 SVG、300 dpi PNG。

## 主要能力

- 遵循 L-α-氨基酸 Fischer 投影的立体化学约定。
- 显式绘制非环重原子及其化学键，并核对双键、芳香环与杂环连接位点。
- 统一 `C`、`H` 的 Times New Roman 字体、字面尺寸和基线。
- 中文名称使用宋体，英文缩写使用 Times New Roman。
- 内置20种氨基酸的默认五组组图版式。
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

## 技能结构

```text
SKILL.md
agents/openai.yaml
references/drawing_rules.md
scripts/draw_amino_acid_fischer.py
```

详细化学和排版规则见 [`references/drawing_rules.md`](references/drawing_rules.md)。
