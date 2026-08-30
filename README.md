# BChemFig：生化教材插图技能

`BChemFig` 是面向中文生物化学与食品生物化学教材的 Codex skill，用于制作和审校科学准确、版式统一、可用于印刷和 Word 排版的生化插图。

## 适用内容

- 氨基酸和其他生化分子的结构式、Fischer 投影式与分组结构图。
- 氨基酸、肽和蛋白质相关反应式。
- 酸碱滴定、氧合、层析校准等数据曲线。
- 等电聚焦、亲和层析、凝胶过滤、透析和超滤等原理或装置图。
- α-螺旋、β-折叠、β-转角和超二级结构模型。
- 基于 PDB 与 PyMOL 的蛋白质三维结构图。
- 序列比较、突变位点和短肽二硫键拓扑图。
- 教学流程图、教材数据表和 Word 可编辑公式。

## 统一规范

- 科学正确性优先于与参考图的视觉相似度。
- 中文使用宋体，正式标题和关键标签采用加粗宋体；英文、数字、化学式和生物学缩写使用 Times New Roman。
- 以黑、白、灰为主，控制标注数量、箭头长度和无信息留白。
- 默认输出 SVG 和 300-dpi PNG，并保留白底与透明背景版本。
- 在原始分辨率下逐项核验结构、价键、数据、单位、标签、字体、重叠和裁切。

## 安装

```powershell
git clone git@github.com:libmelo/BChemFig.git "$env:USERPROFILE\.codex\skills\bchemfig"
```

技能调用名：

```text
$bchemfig
```

## 内置确定性脚本

安装 Python 依赖：

```powershell
python -m pip install -r requirements.txt
```

生成20种标准氨基酸 Fischer 投影式：

```powershell
python scripts/draw_amino_acid_fischer.py --output-dir output
```

生成五种氨基酸典型反应式：

```powershell
python scripts/draw_amino_acid_reactions.py --output-dir output-reactions
```

Windows 默认使用 `C:\Windows\Fonts\times.ttf` 和 `C:\Windows\Fonts\simsun.ttc`。其他系统可通过脚本参数指定兼容字体。

## 技能结构

```text
SKILL.md
agents/openai.yaml
references/
  common_standards.md
  drawing_rules.md
  reaction_scheme_rules.md
  data_plots_tables_equations.md
  apparatus_and_flowcharts.md
  protein_structure_models.md
  sequence_and_topology.md
scripts/
  draw_amino_acid_fischer.py
  draw_amino_acid_reactions.py
```
