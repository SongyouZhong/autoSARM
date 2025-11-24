# autoSARM

**Automatic Structure-Activity Relationship Matrix Generator**

一个用于药物化学研究的自动化工具，专门用于生成和分析构效关系（SAR）表格和分子衍生树。

[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org/)
[![RDKit](https://img.shields.io/badge/RDKit-Required-green.svg)](https://www.rdkit.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## 📋 目录

- [项目简介](#项目简介)
- [核心功能](#核心功能)
- [项目结构](#项目结构)
- [安装依赖](#安装依赖)
- [快速开始](#快速开始)
- [使用说明](#使用说明)
- [应用场景](#应用场景)
- [示例](#示例)
- [技术特点](#技术特点)
- [贡献指南](#贡献指南)
- [许可证](#许可证)

---

## 🎯 项目简介

autoSARM 是一个专业的药物化学计算工具，旨在帮助药物化学家和研究人员快速分析和可视化化合物系列的构效关系。通过自动化的分子片段化和统计分析，该工具能够：

- 🔬 自动识别分子核心骨架和取代基
- 📊 生成一维和二维 SAR 表格
- 🌳 构建分子衍生关系树
- 📈 可视化展示构效关系模式
- 💊 辅助先导化合物优化

该项目特别适用于**激酶抑制剂**、**GPCR 配体**等系列化合物的研发分析。

---

## ✨ 核心功能

### 1. 自动化分子片段化

这是 autoSARM 的核心功能之一，用于将完整分子智能拆分为核心骨架和取代基片段。

#### 1.1 智能切割算法

**基于化学规则的片段化策略**：

- **环旁边的单键切割**：
  - 识别所有环系统（包括芳香环、脂肪环、杂环）
  - 切割与环直接相连的单键
  - 保护环结构完整性
  - 示例：`Ph-CH2-R` → `Ph-*` + `*-CH2-R`

- **环-环之间的单键切割**（可选）：
  - 识别联苯、萘类等多环系统
  - 切割连接两个环的单键
  - 适用于双芳环骨架优化
  - 示例：`Ph-Ph-R` → `Ph-*` + `*-Ph-R`

- **两轮迭代片段化策略**：
  - **第一轮**：识别主要核心骨架和主要取代基
  - **第二轮**：进一步细分大的取代基
  - 自动去重和频率统计
  - 只保留出现频次 ≥ 阈值的片段（避免噪声）

**输入文件要求**：
```csv
smiles,IC50_uM,compound_id
Cc1ccc(C2CCN(C(=O)c3ccccc3)CC2)cc1,0.15,COMP001
Cc1ccc(C2CCN(C(=O)c3cccc(F)c3)CC2)cc1,0.08,COMP002
```
- 必须包含 `smiles` 列
- 活性数据列（如 `IC50_uM`）可以有多个

**参数控制**：
```python
# 在 utils/sarm_utils.py 中的 frag_mol_near_ring 函数
pos_args = {
    'RR': True,      # 环-环之间单键切割
    'nRnR': False    # 非环单键切割（一般不推荐）
}
```

**输出结果**：
- `Frag_round1_count.csv`：第一轮片段及出现次数
- `Frag_round2_count.csv`：第二轮细分片段及出现次数
- 每个片段用带 `*` 的 SMILES 表示（`*` 为连接点）

**使用场景**：
- 大规模化合物库的自动解析
- 识别系列化合物的共同核心结构
- 为 SAR 分析准备片段数据

#### 1.2 核心结构识别

**功能描述**：
- **Murcko scaffold 提取**：自动提取分子的母核骨架
- **药效团识别**：保留关键的杂原子和官能团
- **R 基团自动分类**：
  - 识别所有取代位点
  - 标记取代位置（Left, Right, Center 等）
  - 统计每个位置的取代基种类和频次

**算法细节**：
```python
# 示例：从分子提取核心
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

mol = Chem.MolFromSmiles('Cc1ccc(C2CCN(C(=O)c3ccccc3)CC2)cc1')
scaffold = MurckoScaffold.GetScaffoldForMol(mol)
# 结果：c1ccc(C2CCNCC2)cc1 (哌啶-苯基核心)
```

**使用场景**：
- 专利分析：识别权利要求的核心结构
- 化合物聚类：按骨架类型分组
- 骨架跃迁：探索不同核心结构的活性对比

### 2. SAR 表格生成

自动生成一维和二维构效关系表格，直观展示取代基对活性的影响。

#### 2.1 1D SAR 表格（单点取代分析）

**功能描述**：
分析单个取代位点的构效关系，生成三种类型的表格：

- **Left Table（左侧取代基分析）**：
  - 固定右侧片段，变化左侧取代基
  - 每行代表一个左侧 R 基团
  - 每列代表一个右侧核心片段
  - 单元格显示化合物活性统计

- **Right Table（右侧取代基分析）**：
  - 固定左侧片段，变化右侧取代基
  - 每行代表一个右侧 R 基团
  - 每列代表一个左侧核心片段

- **Single-cut Table（整体分子比对）**：
  - 单切点片段化（只切一个键）
  - 适用于核心骨架变化分析
  - 用于比较结构相似但骨架不同的化合物

**输入文件**：
```csv
smiles,IC50_JAK1,IC50_JAK3,Selectivity
# 系列化合物，共享相似核心骨架
```

**输出文件**：
```
Left_Table/
├── Table_0_Left.csv      # 第一个核心的左侧取代表
├── Table_1_Left.csv
└── ...

Right_Table/
├── Table_0_Right.csv     # 第一个核心的右侧取代表
├── Table_1_Right.csv
└── ...
```

**表格内容**：
```
        | Core_1        | Core_2        | Core_3
--------|---------------|---------------|---------------
R1: -CH3|  IC50=12.5±3  |  IC50=8.3±2   |  IC50=15.2±4
        |  n=3          |  n=2          |  n=4
--------|---------------|---------------|---------------
R2: -CF3|  IC50=5.6±1   |  IC50=2.9±0.5 |  IC50=6.8±1.2
        |  n=5          |  n=3          |  n=2
```

**使用场景**：
- **单点优化**：优化单个取代位点的 R 基团
- **趋势分析**：观察电子效应、立体效应的影响
- **快速筛选**：识别最优 R 基团候选

#### 2.2 2D SAR 表格（双点取代矩阵）

**功能描述**：
同时分析两个取代位点的构效关系，展示加和效应或协同效应。

- **Combine Table（完整的二维矩阵）**：
  - 行：左侧 R 基团
  - 列：右侧 R 基团
  - 单元格：对应双取代化合物的活性数据
  - 自动合并左右表格
  - 智能去重和排序

**输入文件**：
- 来自第一轮片段化的结果
- 必须有足够的双取代化合物（建议 ≥ 3 个样本/单元格）

**输出文件**：
```
Combine_Table/
├── Table_0_combine.csv    # 第一个核心的组合表
├── Table_1_combine.csv
└── ...
```

**表格结构示例**：
```
              | R_right: -H  | R_right: -F  | R_right: -Cl
--------------|--------------|--------------|-------------
R_left: -CH3  |   12.5 (3)   |    8.3 (2)   |   15.2 (4)
R_left: -CH2CH3|   10.2 (2)   |    6.1 (3)   |   12.8 (2)
R_left: -CF3  |    5.6 (5)   |    2.9 (4)   |    6.8 (3)
```
- 数值：活性均值
- 括号内：样本数量

**矩阵解读**：
- **对角线趋势**：评估双取代的协同效应
- **行/列比较**：固定一侧，观察另一侧的影响
- **活性悬崖**：相邻单元格活性差异显著
- **活性热点**：高活性区域集中的位置

**使用场景**：
- **双点优化**：同时优化两个取代位点
- **协同效应**：发现增效或拮抗的取代组合
- **专利布局**：覆盖活性组合的最优空间

#### 2.3 统计信息自动计算

每个表格单元格包含丰富的统计信息：

- **均值（Mean）**：该组合的平均活性
- **标准差（Std）**：数据离散程度
- **中位数（Median）**：对异常值更稳健
- **最小值/最大值（Min/Max）**：活性范围
- **样本数量（Count）**：数据可靠性指标

**输出格式**（CSV）：
```csv
Left_Fragment,Right_Fragment,Mean_IC50,Std_IC50,Median_IC50,Min_IC50,Max_IC50,Count
*-CH3,Ph-*,12.5,3.2,11.8,8.3,18.2,5
*-CF3,Ph-*,5.6,1.1,5.4,4.2,7.3,8
```

**Excel 增强功能**（使用 `--csv2excel 1`）：
- 自动插入分子结构图
- 条件格式化（活性高低颜色标记）
- 自动排序（按活性或样本数）
- 多 sheet 组织（每个核心一个 sheet）

**使用场景**：
- **数据可靠性评估**：通过 Std 和 Count 判断
- **异常值识别**：Min/Max 远离 Mean
- **优先级排序**：按 Mean 或 Median 排序

#### 2.4 完整生成流程

```bash
# 完整命令示例
python create_sarm.py \
    --csvFile compounds.csv \          # 输入文件
    --column IC50_JAK1 IC50_JAK3 \    # 多个活性列
    --type smiles \                    # 输入类型
    --log 1 \                          # 对活性取对数
    --minimumSite1 3 \                 # 第一轮最小出现次数
    --minimumSite2 3 \                 # 第二轮最小出现次数
    --n_jobs 8 \                       # 并行核心数
    --save_folder SAR_Results \        # 输出文件夹
    --csv2excel 1                      # 生成Excel
```

**处理流程**：
```
1. 读取 CSV 文件并验证
   ↓
2. 第一轮片段化 (环旁边单键)
   ↓
3. 统计片段频次，过滤低频片段
   ↓
4. 第二轮片段化 (可选)
   ↓
5. 生成 Left/Right 单表
   ↓
6. 合并生成 Combine 双表
   ↓
7. 计算统计量并排序
   ↓
8. 输出 CSV 和 Excel (可选)
```

**输出文件总览**：
```
SAR_Results/
├── input.csv                      # 原始输入副本
├── Frag_round1_count.csv          # 第一轮片段统计
├── Frag_round2_count.csv          # 第二轮片段统计
├── Left_Table_info.csv            # 左表索引
├── Right_Table_info.csv           # 右表索引
├── Combine_Table_info.csv         # 组合表索引
├── singleCut_Table_info.csv       # 单切表索引
├── Left_Table/                    # 左侧取代表
├── Right_Table/                   # 右侧取代表
├── Combine_Table/                 # 二维组合表
└── GENERATION_SUMMARY.md          # 生成报告
```

### 3. SAR 树可视化

将分子衍生关系以树状图形式展示，直观呈现化合物优化路径和活性趋势。

#### 3.1 层级结构展示

**功能描述**：
- **基于相似性构建树**：
  - 使用 Morgan 指纹（ECFP4）计算分子相似度
  - Tanimoto 相似度作为距离度量
  - 自动聚类并构建父子关系
  
- **可自定义最大层级**：
  - 通过 `--maxLevel` 参数控制树的深度
  - 默认 5 层（根节点 → 第4代衍生物）
  - 避免过深导致可读性下降

- **支持双切割和单切割片段**：
  - `double-cut`：基于 Combine_Table 的双取代片段
  - `single-cut`：基于 singleCut_Table 的单切代片段
  - 可同时展示两种片段在树中的位置

**输入文件要求**：
```
1. SAR 表格文件（由 create_sarm.py 生成）：
   - Combine_Table_info.csv
   - singleCut_Table_info.csv
   
2. 核心片段 SMILES（带 * 标记连接点）：
   例如：*CN1CCC(c2ccc3[nH]c(-c4cc(CO*)c5ncnn5c4)c(C(C)C)c3c2)CC1
```

**命令示例**：
```bash
python create_tree.py \
    --fragment_core "*CN1CCC(c2ccc3[nH]c(-c4cc(CO*)c5ncnn5c4)c(C(C)C)c3c2)CC1" \
    --rootTitle "Table_100_combine" \
    --workFolder ./SAR_Results \
    --maxLevel 5 \
    --treeContent "['double-cut','single-cut']"
```

**树的构建逻辑**：
```
根节点 (核心骨架)
├── 第1层: 直接衍生物 (单个R基团变化)
│   ├── 化合物A (IC50=10.5 μM)
│   ├── 化合物B (IC50=3.2 μM) ← 活性提升
│   └── 化合物C (IC50=15.8 μM)
│
├── 第2层: 双重修饰衍生物
│   ├── 基于化合物B的衍生
│   │   ├── 化合物D (IC50=1.5 μM) ← 继续优化
│   │   └── 化合物E (IC50=4.1 μM)
│   └── ...
│
└── 第3-5层: 更深层次的优化
```

#### 3.2 智能高亮

**功能描述**：
根据活性数据自动高亮符合条件的化合物，快速识别关键分子。

**高亮条件语法**：
```python
# 通过 --highlightDict 参数指定（JSON 格式）
highlightDict = [
    {
        'col': 'IC50_uM',      # 活性列名
        'type': 'means',        # 统计类型: means, median, min, max
        'relation': '<',        # 关系: <, >, <=, >=, ==
        'value': 1.0            # 阈值
    },
    {
        'col': 'Selectivity',
        'type': 'means',
        'relation': '>',
        'value': 10.0
    }
]
```

**高亮示例**：
```bash
python create_tree.py \
    --fragment_core "核心SMILES" \
    --rootTitle "Table_100_combine" \
    --workFolder ./SAR_Results \
    --highlightDict "[{'col':'IC50_uM', 'type':'means', 'relation':'<', 'value':1.0}]"
```

**视觉效果**：
- ✅ **绿色高亮**：满足所有高亮条件的化合物
- ⚪ **白色**：普通化合物
- 🔵 **蓝色边框**：根节点
- 数字标签：活性值显示在节点旁边

**使用场景**：
- **快速定位先导化合物**：高亮 IC50 < 1 μM
- **选择性筛选**：高亮选择性 > 10 倍
- **多参数优化**：同时满足活性、选择性、溶解度等

#### 3.3 支持多靶点同时展示

**功能描述**：
在树中同时显示多个靶点的活性数据，辅助多靶点药物设计。

**输入数据格式**：
```csv
smiles,IC50_JAK1,IC50_JAK3,Selectivity_JAK3_JAK1
CCO,12.5,125,10.0
c1ccccc1,8.3,150,18.1
```

**树节点显示**：
```
[化合物 ID]
JAK1: 12.5 μM
JAK3: 125 μM
Selectivity: 10x
```

**高亮策略**：
```python
# 双靶点活性优化
highlightDict = [
    {'col': 'IC50_JAK1', 'type': 'means', 'relation': '<', 'value': 10.0},
    {'col': 'IC50_JAK3', 'type': 'means', 'relation': '<', 'value': 100.0}
]

# 选择性优化
highlightDict = [
    {'col': 'Selectivity_JAK3_JAK1', 'type': 'means', 'relation': '>', 'value': 15.0}
]
```

**使用场景**：
- **双靶点抑制剂**：JAK1/JAK3, CDK2/CDK4, TLR7/TLR8
- **选择性优化**：提高靶点选择性，降低脱靶效应
- **ADMET 多参数**：同时考虑活性、溶解度、代谢稳定性

#### 3.4 自定义筛选条件

**功能描述**：
灵活定义筛选条件，只展示感兴趣的化合物子集。

**筛选参数**：
```python
# 在 tree_utils.py 中的 filter_compounds 函数
filter_conditions = {
    'IC50_min': 0.1,           # 最小活性阈值
    'IC50_max': 10.0,          # 最大活性阈值
    'Count_min': 3,            # 最小样本数
    'MolWeight_max': 500       # 分子量限制
}
```

**应用场景**：
- **类药性过滤**：MW < 500, LogP < 5
- **数据质量控制**：只显示 Count ≥ 3 的高可信数据
- **活性范围聚焦**：只看 IC50 在 0.1-10 μM 的化合物

#### 3.5 输出格式

**PDF 树状图**：
- 高质量矢量图（可无限缩放）
- 适合论文发表和演示
- 文件名：`{rootTitle}_tree.pdf`

**PNG 分子结构图**：
- 树中每个节点的分子结构
- 2D 结构式，带原子编号（可选）
- 保存在 `{workFolder}/Trees/structures/`

**文本格式的树结构**：
- 易于阅读的文本树
- 包含所有活性数据
- 文件名：`{rootTitle}_tree.txt`

**示例输出**：
```
根节点: 核心骨架 (Table_100_combine)
│
├── [Level 1] 化合物_001 (IC50: 12.5 μM, n=3)
│   ├── [Level 2] 化合物_045 (IC50: 8.3 μM, n=2) ✓
│   │   └── [Level 3] 化合物_123 (IC50: 1.5 μM, n=4) ✓✓
│   └── [Level 2] 化合物_067 (IC50: 15.2 μM, n=5)
│
└── [Level 1] 化合物_002 (IC50: 18.1 μM, n=2)
    └── [Level 2] 化合物_089 (IC50: 6.7 μM, n=3)

✓ = 满足高亮条件
```

#### 3.6 完整使用流程

```bash
# 步骤 1: 生成 SAR 表格（如果还没有）
python create_sarm.py --csvFile compounds.csv --column IC50_uM --save_folder SAR_Results

# 步骤 2: 选择核心片段（从 Combine_Table_info.csv 中查看）
# 假设选择 Table_100_combine 作为树的根节点

# 步骤 3: 生成 SAR 树
python create_tree.py \
    --fragment_core "*CN1CCC(c2ccc3[nH]c(-c4cc(CO*)c5ncnn5c4)c(C(C)C)c3c2)CC1" \
    --rootTitle "Table_100_combine" \
    --workFolder ./SAR_Results \
    --maxLevel 5 \
    --treeContent "['double-cut','single-cut']" \
    --highlightDict "[{'col':'IC50_uM', 'type':'means', 'relation':'<', 'value':1.0}]"

# 步骤 4: 查看结果
# SAR_Results/Trees/Table_100_combine_tree.pdf
# SAR_Results/Trees/Table_100_combine_tree.txt
```

**使用场景总结**：
- **优化路径可视化**：展示从先导化合物到候选药物的演化过程
- **专利分析**：理解竞争对手的化合物优化策略
- **团队沟通**：向非专业人员解释SAR关系
- **历史数据回顾**：梳理项目中所有合成的化合物

### 4. 基于核心结构的 SAR 分析

针对已知核心骨架，系统化分析所有可能的取代位点和R基团组合。

#### 4.1 单核心分析功能

**功能描述**：
- 指定一个核心结构（core scaffold）
- 自动识别核心上的所有虚拟原子（dummy atoms，用 `*` 表示）
- 枚举所有在该核心上进行取代的化合物
- 生成针对该核心的专属 SAR 表格

**输入文件要求**：
```csv
# CSV 文件（与 create_sarm.py 相同）
smiles,IC50_JAK1,IC50_JAK3,Selectivity
c1ccc2[nH]c(-c3cc(CO)c4ncnn4c3)c(C(C)C)c2c1,5.2,52,10.0
Cc1ccc2[nH]c(-c3cc(CO)c4ncnn4c3)c(C(C)C)c2c1,3.8,68,17.9
...
```

**核心结构定义**：
```python
# 核心 SMILES（不带虚拟原子）
core = "c1ccc2[nH]c(-c3cc(CO)c4ncnn4c3)c(C(C)C)c2c1"

# 或者带虚拟原子标记取代位置
core = "*c1ccc2[nH]c(-c3cc(CO*)c4ncnn4c3)c(C(C)C)c2c1"
```

**命令示例**：
```bash
python print_sar_single_core.py \
    --csvFile compounds.csv \
    --core "c1ccc2[nH]c(-c3cc(CO)c4ncnn4c3)c(C(C)C)c2c1" \
    --actCols IC50_JAK1 IC50_JAK3 Selectivity \
    --n_jobs 8 \
    --save_folder Core_SAR_Results
```

#### 4.2 R 基团枚举

**功能描述**：
- 自动检测核心结构上的所有取代位点
- 提取每个位点的所有 R 基团变化
- 统计每个 R 基团的出现频次
- 按活性或频次排序

**枚举逻辑**：
```python
# 假设核心有 3 个取代位点：R1, R2, R3
核心: *-Core-*-*
       ↓   ↓  ↓
      R1  R2  R3

# 自动提取：
R1_groups = ['-CH3', '-CF3', '-H', '-OCH3', ...]
R2_groups = ['-Ph', '-Py', '-Thienyl', ...]
R3_groups = ['-COOH', '-CONH2', '-CN', ...]
```

**输出文件**：
```
Core_SAR_Results/
├── core_info.csv                 # 核心结构信息
├── R1_groups.csv                 # R1 位置的所有基团
├── R2_groups.csv                 # R2 位置的所有基团
├── R3_groups.csv                 # R3 位置的所有基团
└── full_enumeration.csv          # 完整的 R 基团组合
```

**R 基团统计表格式**：
```csv
R_group,SMILES,Count,Mean_IC50_JAK1,Std_IC50_JAK1,Mean_Selectivity
-CH3,*C,15,8.5,3.2,12.3
-CF3,*C(F)(F)F,8,3.2,1.1,18.7
-OCH3,*OC,12,12.1,4.5,8.9
```

**使用场景**：
- **系统化SAR研究**：覆盖所有可能的取代位点
- **R基团优选**：快速识别最佳取代基
- **专利权利要求撰写**：列举所有活性取代基

#### 4.3 自动提取所有取代位点

**功能描述**：
- 智能识别核心上的取代位置（包括芳香环、杂环、脂肪链）
- 自动命名取代位点（R1, R2, R3, ...）
- 处理多个虚拟原子的复杂核心

**取代位点检测算法**：
```python
from rdkit import Chem

mol = Chem.MolFromSmiles("*c1ccc(*)cc1*")  # 3个取代位点
dummy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
n_sites = len(dummy_atoms)  # 3 个取代位点
```

**复杂核心示例**：
```python
# 核心：吲哚-吡唑并嘧啶双环系统（3个取代位点）
core = "*c1ccc2[nH]c(-c3cc(CO*)c4ncnn4c3)c(C(C)C)c2c1"
#       ↑R1                    ↑R2              (R3在C(C)C位置)

# 自动识别：
# R1: 吲哚5位（芳香环取代）
# R2: 吡唑并嘧啶侧链（羟甲基邻位）
# R3: 吲哚3位（异丙基位置，潜在取代）
```

**使用场景**：
- **多取代位点优化**：同时优化2-4个位点
- **正交设计**：设计组合化学库
- **构效关系模式发现**：识别位点间的相互作用

#### 4.4 带活性数据的 Excel 报告

**功能描述**：
- 自动生成带分子结构图的 Excel 文件
- 每个取代位点一个 sheet
- 包含所有活性数据和统计量
- 支持条件格式化（高活性绿色，低活性红色）

**Excel 文件结构**：
```
Core_SAR_Report.xlsx
├── [Sheet 1] Summary - 总览
│   ├── 核心结构图
│   ├── 取代位点示意图
│   └── 统计摘要
│
├── [Sheet 2] R1_Analysis - R1 基团分析
│   ├── 列A: R基团结构图
│   ├── 列B: R基团 SMILES
│   ├── 列C-F: 活性数据 (IC50_JAK1, IC50_JAK3, Selectivity, ...)
│   └── 列G: 样本数
│
├── [Sheet 3] R2_Analysis - R2 基团分析
├── [Sheet 4] R3_Analysis - R3 基团分析
│
└── [Sheet 5] Full_Matrix - 完整组合矩阵
    └── R1 × R2 × R3 的三维矩阵（分页显示）
```

**自动生成分子结构图**：
```python
from rdkit import Chem
from rdkit.Chem import Draw

# 为每个 R 基团生成结构图
for r_group in r_groups:
    mol = Chem.MolFromSmiles(r_group)
    img = Draw.MolToImage(mol, size=(200, 200))
    # 插入到 Excel 单元格
```

**条件格式化规则**：
```python
# 活性高（绿色）：IC50 < 5 μM
# 活性中等（黄色）：5 μM ≤ IC50 < 20 μM
# 活性低（红色）：IC50 ≥ 20 μM
```

**生成命令**：
```bash
python print_sar_single_core.py \
    --csvFile compounds.csv \
    --core "c1ccc2[nH]c(...)c2c1" \
    --actCols IC50_JAK1 IC50_JAK3 Selectivity \
    --output_excel 1 \              # 启用 Excel 输出
    --excel_name "JAK_Core_SAR.xlsx"
```

**使用场景**：
- **项目汇报**：生成专业的 SAR 分析报告
- **数据共享**：与团队成员分享结构化数据
- **决策支持**：快速浏览和比较不同 R 基团

#### 4.5 完整工作流程

```bash
# 步骤 1: 准备包含核心骨架系列的化合物数据
# compounds.csv 中应包含多个共享相同核心的化合物

# 步骤 2: 确定核心结构
# 方法1：手动定义（已知核心）
core="c1ccc2[nH]c(-c3cc(CO)c4ncnn4c3)c(C(C)C)c2c1"

# 方法2：自动提取（使用 Murcko scaffold）
python -c "
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
mol = Chem.MolFromSmiles('your_example_compound')
core = MurckoScaffold.GetScaffoldForMol(mol)
print(Chem.MolToSmiles(core))
"

# 步骤 3: 运行单核心 SAR 分析
python print_sar_single_core.py \
    --csvFile compounds.csv \
    --core "$core" \
    --actCols IC50_JAK1 IC50_JAK3 Selectivity \
    --n_jobs 8 \
    --save_folder JAK_Core_SAR \
    --output_excel 1

# 步骤 4: 查看结果
# JAK_Core_SAR/Core_SAR_Report.xlsx
# JAK_Core_SAR/R1_groups.csv
# JAK_Core_SAR/R2_groups.csv
# ...
```

**实际应用案例**：

**案例1：激酶抑制剂优化**
```python
# 核心：ATP 竞争性激酶抑制剂骨架
core = "c1cnc2ncnc(N)c2c1"  # 嘌呤核心

# 分析：
# - R1 (6位)：铰链区氢键供体/受体
# - R2 (9位)：疏水口袋取代
# - R3 (2位)：选择性口袋探索
```

**案例2：GPCR 配体设计**
```python
# 核心：哌嗪-苯并咪唑
core = "c1ccc2nc([nH]c2c1)N3CCNCC3"

# 分析：
# - R1 (苯并咪唑2位)：与受体 TM3 相互作用
# - R2 (哌嗪N端)：溶剂暴露区域，影响选择性
```

**使用场景总结**：
- **聚焦优化**：针对特定核心的深度优化
- **专利设计**：系统化列举核心的所有变体
- **文献调研**：提取文献中某一核心的所有化合物
- **虚拟筛选**：基于核心的虚拟枚举和打分

### 5. 3D 空间位置分析（基于结构的片段优化）

这是一个强大的功能模块，专门用于分析分子片段在蛋白质结合口袋中的空间分布，辅助基于结构的药物设计。

#### 5.1 对接结果处理

**功能描述**：
- 从 Schrödinger Glide 对接结果（maegz 格式）批量提取 SDF 文件
- 支持按对接评分（GlideScore）自动筛选高质量构象
- 批量处理多个对接姿态，保留最优结合模式

**输入文件要求**：
```
文件格式：MAEGZ (Maestro compressed format)
来源：Schrödinger Glide 或其他 Schrödinger 对接工具
依赖：需要安装 Schrödinger Suite (2021-2 或更高版本)
```

**使用场景**：
- 处理虚拟筛选后的大量对接结果
- 提取高评分化合物的 3D 构象用于进一步分析
- 为网格化分析准备标准化的 SDF 输入文件

**命令示例**：
```bash
python scripts/extract_sdf_fromMaegz.py \
    --maegz docking_results.maegz \
    --saveDir ./sdf_output \
    --scoreLimit 6.0 \
    --format sd
```

#### 5.2 原子位置网格化分析

**功能描述**：
- **3D 空间网格划分**：将结合口袋区域划分为规则的 3D 网格单元（默认 0.5 Å × 0.5 Å × 0.5 Å）
- **原子密度统计**：统计每个网格单元中不同类型原子的出现频率
- **原子类型分类**：
  - 按元素类型：C（碳）、N（氮）、O（氧）、S（硫）、F/Cl/Br（卤素）
  - 按芳香性：芳香碳 vs 非芳香碳
  - 按环状态：环状原子特殊标记
- **热点识别**：自动识别高密度区域（片段优先占据位置）

**输入文件要求**：
```
文件格式：SDF (Structure Data File)
内容要求：必须包含 3D 坐标信息
推荐来源：分子对接软件输出（AutoDock Vina, Glide, GOLD, MOE 等）
文件组织：所有 SDF 文件放在同一文件夹中
```

**分析流程**：
```
1. 读取 SDF 文件夹中的所有分子构象
   ↓
2. 提取每个原子的 3D 坐标 (x, y, z)
   ↓
3. 根据坐标将原子分配到对应的网格单元
   ↓
4. 统计每个网格单元的原子类型和数量
   ↓
5. 输出网格密度统计 CSV 文件
   ↓
6. （可选）生成 PyMOL 可视化脚本
```

**输出文件**：
- `atom_count.csv`：包含每个网格单元的原子统计信息
  - 列：`grid_x, grid_y, grid_z, C_count, N_count, O_count, aromatic_C, aliphatic_C, ...`
- `grid_visualization.pse`：PyMOL 会话文件（可视化）

**使用场景**：
- **热点分析**：识别结合口袋中配体原子高频出现的区域
- **片段设计指导**：了解哪些位置适合放置特定类型的原子或基团
- **构效关系解释**：理解为什么某些取代基活性更好（空间匹配度）
- **虚拟筛选验证**：验证对接姿态的合理性

#### 5.3 片段空间位置评分

**功能描述**：
- **片段空间占据计算**：量化片段在结合口袋中占据的 3D 空间体积
- **相对位置评分**：
  - 与参考片段（如已知抑制剂）的空间重叠度
  - 片段原子在高密度网格区域的分布比例
  - 偏离核心结合区域的惩罚评分
- **片段优先级排序**：根据空间匹配度自动排序候选片段

**输入文件要求**：
```
1. 片段 SDF 文件：待评分的片段 3D 构象
2. 参考 SDF 文件（可选）：已知活性化合物的结合构象
3. 网格密度文件：由网格化分析生成的 atom_count.csv
```

**评分算法**：
```python
# 伪代码示例
score = (
    0.4 × 高密度区域占据率 +
    0.3 × 与参考片段的空间重叠度 +
    0.2 × 关键相互作用位点匹配度 +
    0.1 × 片段大小适配度
)
```

**使用场景**：
- **片段筛选**：从大量片段中筛选出最适合结合口袋的候选
- **R 基团优化**：比较不同 R 基团在特定取代位置的空间适配性
- **骨架跃迁**：评估新骨架与已知活性化合物的空间相似性
- **先导优化**：量化评估结构修饰后的空间互补性变化

#### 5.4 PyMOL 3D 可视化

**功能描述**：
- **原子密度热图**：用球体或立方体表示网格单元，颜色深浅表示原子密度
- **多颜色方案**：
  - 按原子类型着色：C (青色), N (蓝色), O (红色), S (黄色)
  - 按芳香性着色：芳香碳 (绿色) vs 脂肪碳 (灰色)
  - 按密度梯度：低密度 (透明) → 高密度 (不透明)
- **叠加显示**：
  - 蛋白质结合口袋（cartoon 或 surface）
  - 已知配体（stick）
  - 原子密度网格（sphere）
  - 候选片段（ball-and-stick）

**生成方法**：
```python
from utils.grid_pos_utils import create_grid_pse

create_grid_pse(
    csvFile='atom_count.csv',      # 网格密度文件
    outputPse='visualization.pse',  # 输出 PyMOL 文件
    group='atomDensity',            # 对象分组名称
    colorScheme='by_atom_type'      # 颜色方案
)
```

**可视化示例**：
```
PyMOL 中加载后显示：
- 蛋白质表面（半透明灰色）
- 结合口袋网格热图（彩色球体）
  - 红色高密度区域 = 配体氧原子常出现位置
  - 蓝色高密度区域 = 配体氮原子常出现位置
  - 青色高密度区域 = 配体碳骨架区域
- 参考配体（绿色 stick）
```

**使用场景**：
- **结果展示**：为论文、报告、演讲制作高质量 3D 图
- **交互式探索**：在 PyMOL 中旋转、缩放，从不同角度观察
- **假说验证**：直观验证片段设计的合理性
- **团队讨论**：与生物学家、结构生物学家讨论设计策略

#### 5.5 完整工作流程示例

```python
# 步骤 1: 从对接结果提取 SDF
from utils.grid_pos_utils import extract_sdf
extract_sdf('docking.maegz', './sdf_output', scoreLimit=6.0)

# 步骤 2: 网格化分析
df_atmCount, grid_params = grid_atom('./sdf_output')
df_atmCount.to_csv('atom_density.csv', index=False)

# 步骤 3: 计算片段评分
from utils.grid_pos_utils import score_fragment_position
scores = score_fragment_position(
    fragment_sdf='my_fragment.sdf',
    grid_csv='atom_density.csv',
    reference_sdf='known_inhibitor.sdf'
)

# 步骤 4: 生成 PyMOL 可视化
create_grid_pse('atom_density.csv', 'final_visualization.pse')
```

**实际应用案例**：
- **激酶抑制剂设计**：分析铰链区氢键供体/受体的最优位置
- **GPCR 配体优化**：识别疏水口袋和极性相互作用热点
- **蛋白酶抑制剂**：优化 P1、P2、P3 位置的侧链取代基

---

## 📁 项目结构

```
autoSARM/
├── create_sarm.py              # 主程序：创建 SAR 表格
├── create_tree.py              # 主程序：生成 SAR 树
├── print_sar_single_core.py    # 基于单核心结构生成 SAR（支持多个 R 基团）
├── utils/                      # 核心工具模块
│   ├── __init__.py
│   ├── sarm_utils.py          # SAR 核心算法（片段化、表格生成）
│   ├── tree_utils.py          # 树结构构建和可视化
│   ├── dual_target_utils.py   # 双靶点分析和片段匹配
│   ├── grid_pos_utils.py      # 3D 空间位置分析和网格化
│   └── common_utils.py        # 通用工具函数（分子处理、Excel 导出）
├── scripts/                    # 辅助脚本
│   ├── extract_sdf_fromMaegz.py    # 从 Schrödinger 对接结果提取 SDF
│   └── BindingDB_processing.ipynb   # BindingDB 数据处理
├── CDK1-2/                     # 案例研究：CDK 抑制剂
│   ├── Create_tree_AIXB7.ipynb
│   ├── fragmentize_CDK2.ipynb
│   └── CDK2_SMILES/           # 包含完整的片段化结果和 SAR 表格
├── create_tree.ipynb           # Jupyter 示例：TLR7/8 树状图生成
├── print_sar.ipynb             # Jupyter 示例：SAR 表格生成
└── test_functions.ipynb        # 功能测试
```

---

## 🔧 安装依赖

### 系统要求

- Python 3.7+
- Linux/macOS/Windows

### 核心依赖

```bash
# 化学信息学核心库
rdkit>=2020.09
```

### 数据处理和可视化

```bash
# 数据处理
pandas>=1.0.0
numpy>=1.18.0

# 并行计算
pandarallel>=1.5.0

# 可视化
matplotlib>=3.0.0
seaborn>=0.10.0
plotly>=4.0.0

# Excel 支持
openpyxl>=3.0.0

# 图形生成
graphviz>=0.14.0  # 需要系统级 Graphviz
svglib>=1.0.0
reportlab>=3.5.0
IPython
pip install graphviz
micromamba install graphviz
# 分子可视化（用于 3D 空间分析）
pymol-open-source>=2.5.0  # 可选，用于 3D 网格可视化
```

### Schrödinger 套件（可选）

用于对接结果处理：
```bash
# 需要安装 Schrödinger Suite（2021-2 或更高版本）
# 用于从 maegz 文件提取 SDF 结构
```

### 安装步骤

1. **克隆仓库**

```bash
git clone git@github.com:SongyouZhong/autoSARM.git
cd autoSARM
```

2. **创建虚拟环境**（推荐）

```bash
conda create -n autosarm python=3.8
conda activate autosarm
```

3. **安装 RDKit**

```bash
conda install -c conda-forge rdkit
```

4. **安装其他依赖**

```bash
pip install pandas numpy pandarallel matplotlib seaborn plotly openpyxl svglib reportlab
```

5. **安装 Graphviz**（用于树状图生成）

```bash
# Ubuntu/Debian
sudo apt-get install graphviz

# macOS
brew install graphviz

# Windows
# 从 https://graphviz.org/download/ 下载安装
```

```bash
pip install graphviz
```

---

## 🚀 快速开始

### 示例 1：生成 SAR 表格

准备输入文件 `compounds.csv`：

```csv
smiles,IC50_uM,Selectivity
Cc1ccc(C2CCN(C(=O)c3ccccc3)CC2)cc1,0.15,10.5
Cc1ccc(C2CCN(C(=O)c3cccc(F)c3)CC2)cc1,0.08,15.2
...
```

运行命令：

```bash
python create_sarm.py \
    --csvFile compounds.csv \
    --column IC50_uM Selectivity \
    --type smiles \
    --log 1 \
    --minimumSite1 3 \
    --minimumSite2 3 \
    --n_jobs 8 \
    --save_folder SAR_Results \
    --csv2excel 1

python create_sarm.py \
  --csvFile SAR_Results/input.csv \
  --column IC50_uM \
  --type smiles \
  --log 1 \
  --minimumSite1 3 \
  --minimumSite2 3 \
  --n_jobs 8 \
  --save_folder SAR_Results_new \
  --csv2excel 1
```

**输出文件**：
- `SAR_Results/Combine_Table_info.csv` - 主表信息
- `SAR_Results/Combine_Table/` - 完整的 2D SAR 表格
- `SAR_Results/Left_Table/` - 左侧取代基表格
- `SAR_Results/Right_Table/` - 右侧取代基表格

### 示例 2：生成 SAR 树

```bash
python create_tree.py \
    --fragment_core "*CN1CCC(c2ccc3[nH]c(-c4cc(CO*)c5ncnn5c4)c(C(C)C)c3c2)CC1" \
    --rootTitle "Table_100_combine" \
    --workFolder ./SAR_Results \
    --maxLevel 5 \
    --treeContent "['double-cut','single-cut']" \
    --highlightDict "[{'col':'IC50_uM', 'type':'means', 'relation':'<', 'value':1.0}]"
```

**输出**：
- PDF 树状图
- 分子结构图片
- 文本格式的树结构

### 示例 3：基于核心结构的 SAR 分析

```bash
python print_sar_single_core.py \
    --csvFile compounds.csv \
    --core "c1ccc2[nH]c(-c3cc(CO)c4ncnn4c3)c(C(C)C)c2c1" \
    --actCols IC50_uM Selectivity \
    --n_jobs 8 \
    --save_folder Core_SAR_Results
```

**输出**：
- 核心结构信息表
- 所有 R 基团枚举
- 带活性数据的 Excel 报告

### 示例 4：3D 对接结果分析（可选功能）

这是一个高级功能，用于基于结构的药物设计。

#### 4.1 从 Schrödinger 对接结果提取 SDF

**前提条件**：
- 已安装 Schrödinger Suite (2021-2 或更高版本)
- 配置环境变量：`SCHRODINGER=/path/to/schrodinger`

**输入文件**：
```
docking_results.maegz  - Glide 对接输出文件
包含内容：
  - 对接的化合物姿态
  - GlideScore 评分
  - 蛋白-配体相互作用信息
```

**提取命令**：
```bash
python scripts/extract_sdf_fromMaegz.py \
    --maegz docking_results.maegz \
    --saveDir ./sdf_output \
    --scoreLimit 6.0 \
    --format sd
```

**参数说明**：
- `--maegz`: Schrödinger 对接输出文件（必需）
- `--scoreLimit`: 最小对接评分阈值（默认 6.0）
  - GlideScore > 6.0 通常表示较好的结合
  - 可根据具体靶点调整（推荐范围 5.0-8.0）
- `--format`: 输出格式
  - `sd`: SDF 格式（推荐，通用性好）
  - `mol2`: Mol2 格式（包含原子类型信息）
- `--saveDir`: SDF 文件保存目录

**输出结果**：
```
sdf_output/
├── compound_001_pose1.sdf
├── compound_002_pose1.sdf
├── compound_003_pose1.sdf
└── ...
```

**使用场景**：
- 虚拟筛选后处理大量对接姿态
- 为网格化分析准备标准输入文件
- 提取高评分化合物用于进一步研究

#### 4.2 原子位置网格化分析

**功能**：将结合口袋划分为 3D 网格，统计每个网格单元的原子密度。

**Python API 示例**：
```python
from utils.grid_pos_utils import grid_atom

# 网格化分析
df_atmCount, grid_params = grid_atom(
    sdfFolder='./sdf_output',      # SDF 文件夹
    gridSize=0.5,                   # 网格单元大小 (Å)
    atomTypes=['C', 'N', 'O', 'S', 'F', 'Cl', 'Br']  # 要统计的原子类型
)

# 保存结果
df_atmCount.to_csv('atom_density.csv', index=False)
print(f"网格参数: {grid_params}")
```

**输出文件**：`atom_density.csv`
```csv
grid_x,grid_y,grid_z,C_count,N_count,O_count,aromatic_C,aliphatic_C,ring_atoms,total_count
10.0,15.5,22.3,45,8,12,30,15,50,65
10.5,15.5,22.3,38,5,10,25,13,42,53
...
```

**列含义**：
- `grid_x/y/z`: 网格单元的中心坐标 (Å)
- `C_count`: 碳原子数量
- `N_count`: 氮原子数量
- `O_count`: 氧原子数量
- `aromatic_C`: 芳香碳数量
- `aliphatic_C`: 脂肪碳数量
- `ring_atoms`: 环状原子总数
- `total_count`: 该网格单元的原子总数

**使用场景**：
- **热点识别**：找出配体原子高频出现的区域
- **相互作用分析**：
  - 高氮/氧密度区域 → 氢键相互作用热点
  - 高碳密度区域 → 疏水相互作用区域
- **片段设计指导**：了解在特定位置适合放置哪种原子

#### 4.3 片段空间位置评分

**功能**：量化评估片段在结合口袋中的空间匹配度。

**Python API 示例**：
```python
from utils.grid_pos_utils import score_fragment_position

# 计算片段评分
scores = score_fragment_position(
    fragment_sdf='my_fragment.sdf',       # 待评分的片段
    grid_csv='atom_density.csv',          # 网格密度文件
    reference_sdf='known_inhibitor.sdf',  # 参考化合物（可选）
    scoring_method='density_overlap'       # 评分方法
)

print(f"片段空间匹配得分: {scores['total_score']:.2f}")
print(f"高密度区域占据率: {scores['hotspot_coverage']:.1f}%")
print(f"与参考化合物重叠度: {scores['reference_overlap']:.1f}%")
```

**批量评分示例**：
```python
import pandas as pd
from glob import glob

# 批量评分多个片段
results = []
for frag_sdf in glob('fragments/*.sdf'):
    score = score_fragment_position(
        fragment_sdf=frag_sdf,
        grid_csv='atom_density.csv'
    )
    results.append({
        'fragment': frag_sdf,
        'score': score['total_score'],
        'hotspot_coverage': score['hotspot_coverage']
    })

# 保存结果
df_scores = pd.DataFrame(results)
df_scores = df_scores.sort_values('score', ascending=False)
df_scores.to_csv('fragment_scores.csv', index=False)

# 查看Top 10片段
print(df_scores.head(10))
```

**使用场景**：
- **片段筛选**：从库中选择最匹配结合口袋的片段
- **R基团优化**：比较不同R基团的空间适配性
- **设计验证**：在合成前预测新设计的合理性

#### 4.4 生成 PyMOL 可视化

**功能**：生成可在 PyMOL 中交互式查看的 3D 原子密度热图。

**Python API 示例**：
```python
from utils.grid_pos_utils import create_grid_pse

# 生成 PyMOL 会话文件
create_grid_pse(
    csvFile='atom_density.csv',
    outputPse='visualization.pse',
    group='atomDensity',
    colorScheme='by_atom_type',      # 颜色方案
    densityThreshold=5               # 只显示原子数 ≥ 5 的网格
)

print("PyMOL 文件已生成: visualization.pse")
print("使用 PyMOL 打开: pymol visualization.pse")
```

**颜色方案选项**：
- `by_atom_type`: 按原子类型着色（C: 青色, N: 蓝色, O: 红色）
- `by_density`: 按密度梯度着色（低密度透明 → 高密度不透明）
- `by_aromaticity`: 按芳香性着色（芳香碳: 绿色, 脂肪碳: 灰色）

**使用场景**：
- **结果展示**：制作高质量 3D 图用于论文、报告
- **交互式探索**：旋转、缩放，从不同角度观察
- **设计讨论**：与团队成员讨论片段设计策略

#### 4.5 完整工作流程示例

```python
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""完整的3D空间分析工作流程"""

from utils.grid_pos_utils import extract_sdf, grid_atom, score_fragment_position, create_grid_pse
import pandas as pd
from glob import glob

# 步骤 1: 从对接结果提取 SDF
print("步骤 1: 提取SDF文件...")
extract_sdf(
    maegz='kinase_docking_results.maegz',
    saveDir='./sdf_extracted',
    scoreLimit=6.5
)

# 步骤 2: 网格化分析
print("步骤 2: 网格化分析...")
df_atmCount, grid_params = grid_atom(
    sdfFolder='./sdf_extracted',
    gridSize=0.5
)
df_atmCount.to_csv('kinase_atom_density.csv', index=False)

# 步骤 3: 评分候选R基团
print("步骤 3: 评分候选R基团...")
r_group_scores = []
for r_sdf in glob('r_groups/*.sdf'):
    score = score_fragment_position(
        fragment_sdf=r_sdf,
        grid_csv='kinase_atom_density.csv',
        reference_sdf='known_inhibitor.sdf'
    )
    r_group_scores.append({
        'R_group': r_sdf.split('/')[-1].replace('.sdf', ''),
        'score': score['total_score'],
        'hotspot_coverage': score['hotspot_coverage']
    })

df_scores = pd.DataFrame(r_group_scores)
df_scores = df_scores.sort_values('score', ascending=False)
df_scores.to_csv('r_group_ranking.csv', index=False)
print("\n=== Top 5 R基团 ===")
print(df_scores.head())

# 步骤 4: 生成可视化
print("\n步骤 4: 生成PyMOL可视化...")
create_grid_pse(
    csvFile='kinase_atom_density.csv',
    outputPse='kinase_hotspots.pse',
    group='Hotspots',
    colorScheme='by_atom_type',
    densityThreshold=10
)

print("\n✅ 分析完成!")
print("请使用 PyMOL 打开 kinase_hotspots.pse 查看结果")
```

---

## 📖 使用说明

### 参数说明

#### `create_sarm.py` 参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--csvFile` | 输入 CSV 文件路径（必需包含 smiles 列） | 必需 |
| `--type` | 分子类型：`smiles` 或 `scaffold` | `smiles` |
| `--column` | 活性数据列名（可多个） | 必需 |
| `--log` | 是否对活性值取对数（0/1） | 0 |
| `--minimumSite1` | 第一轮片段最小出现次数 | 3 |
| `--minimumSite2` | 第二轮片段最小出现次数 | 3 |
| `--n_jobs` | 并行计算核心数 | 8 |
| `--save_folder` | 输出文件夹 | `SAR_Tables` |
| `--csv2excel` | 是否生成 Excel（0/1） | 0 |

#### `create_tree.py` 参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--fragment_core` | 核心片段 SMILES（带 * 标记） | 必需 |
| `--rootTitle` | 树的根节点标题 | 必需 |
| `--workFolder` | SAR 结果文件夹 | 必需 |
| `--maxLevel` | 最大层级 | 5 |
| `--treeContent` | 树内容类型列表 | `['double-cut']` |
| `--highlightDict` | 高亮条件（JSON 格式） | `''` |

#### `print_sar_single_core.py` 参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--csvFile` | 输入 CSV 文件路径 | 必需 |
| `--core` | 核心结构 SMILES | 必需 |
| `--actCols` | 活性数据列名（可多个） | `['ITK','JAK3','selectivity']` |
| `--n_jobs` | 并行计算核心数 | 8 |
| `--save_folder` | 输出文件夹 | `Core_SAR` |

#### `extract_sdf_fromMaegz.py` 参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--maegz` | Schrödinger 对接结果文件 | 必需 |
| `--saveDir` | SDF 输出文件夹 | 必需 |
| `--scoreLimit` | 最小对接评分 | 6.0 |
| `--format` | 输出格式（`sd` 或 `mol2`） | `sd` |

### 输入数据格式

CSV 文件必须包含以下列：
- `smiles`：SMILES 格式的分子结构
- 活性数据列：如 `IC50_uM`, `Ki_nM`, `Selectivity` 等

示例：
```csv
smiles,IC50_JAK1,IC50_JAK2,Selectivity
Cc1ccc(C2CCN(C(=O)c3ccccc3)CC2)cc1,12.5,125,10.0
Cc1ccc(C2CCN(C(=O)c3cccc(F)c3)CC2)cc1,8.3,150,18.1
```

---

## 🎓 应用场景

### 1. 药物发现与先导化合物优化

#### 1.1 先导化合物优化工作流程

**阶段1：初始命中物评估**
```bash
# 输入：高通量筛选（HTS）或虚拟筛选的命中化合物
python create_sarm.py \
    --csvFile HTS_hits.csv \
    --column IC50_uM \
    --minimumSite1 2 \      # 降低阈值以包含更多初始片段
    --save_folder Initial_SAR
```

**使用目的**：
- 快速识别命中物的共同核心骨架
- 发现初步的构效关系模式
- 筛选出值得进一步优化的系列

**阶段2：先导化合物系列扩展**
```bash
# 输入：基于命中物设计合成的系列化合物
python create_sarm.py \
    --csvFile lead_series.csv \
    --column IC50_uM Selectivity Solubility \
    --minimumSite1 3 \
    --minimumSite2 3 \
    --csv2excel 1
```

**使用目的**：
- 生成详细的 2D SAR 矩阵
- 识别活性悬崖（相似结构但活性差异大）
- 发现增效或拮抗的取代组合

**阶段3：先导化合物精细优化**
```bash
# 针对最优核心进行深度优化
python print_sar_single_core.py \
    --csvFile optimized_compounds.csv \
    --core "最优核心SMILES" \
    --actCols IC50 Selectivity LogD Solubility \
    --output_excel 1
```

**使用目的**：
- 系统化枚举所有 R 基团变化
- 多参数优化（活性、选择性、成药性）
- 生成候选药物报告

#### 1.2 快速识别关键药效团

**输入数据**：
```csv
smiles,IC50_uM,notes
# 包含活性化合物和非活性化合物
CCO,0.5,active
CCC,50.0,inactive
```

**分析方法**：
```bash
# 生成 SAR 表格
python create_sarm.py --csvFile compounds.csv --column IC50_uM

# 比较活性和非活性化合物的片段差异
# 在 Combine_Table 中识别：
# - 只在活性化合物中出现的片段 = 必需药效团
# - 导致活性丧失的片段 = 有害基团
```

**使用目的**：
- 定义最小药效团（minimal pharmacophore）
- 识别禁用基团（liability groups）
- 指导后续化合物设计

#### 1.3 SAR 模式分析

**活性悬崖识别**：
```python
# 在 2D SAR 表格中寻找相邻单元格活性差异 > 10倍
# 示例：
#         | R2: -H    | R2: -F
# --------|-----------|----------
# R1: -CH3|  10.5 μM  | 0.8 μM ← 活性悬崖（13倍差异）
# R1: -CF3|  5.2 μM   | 4.9 μM
```

**使用目的**：
- 识别对活性至关重要的取代位点
- 发现微小结构变化导致的活性巨变
- 避免在不敏感位点浪费合成资源

**活性热点发现**：
```python
# 在 SAR 树中追踪活性提升路径
根节点 (IC50=50 μM)
  → 第1代 (IC50=10 μM, 5倍提升)
    → 第2代 (IC50=2 μM, 5倍提升)  ← 持续优化路径
      → 第3代 (IC50=0.5 μM, 4倍提升)
```

**使用目的**：
- 可视化优化历程
- 理解哪些修饰策略最有效
- 指导下一轮设计方向

#### 1.4 系列化合物评估

**批量专利化合物分析**：
```bash
# 输入：从竞争对手专利中提取的化合物
python create_sarm.py \
    --csvFile patent_compounds.csv \
    --column predicted_IC50 \     # 使用预测活性
    --save_folder Patent_Analysis
```

**使用目的**：
- 理解竞争对手的设计策略
- 识别专利覆盖的核心骨架和 R 基团范围
- 寻找设计空间（design around）

---

### 2. 激酶抑制剂设计（专业应用）

autoSARM 在激酶抑制剂领域具有特别优势，已成功应用于多个激酶靶点。

#### 2.1 JAK 激酶抑制剂（JAK1/JAK2/JAK3）

**研究目标**：
- 开发高选择性 JAK3 抑制剂（治疗类风湿性关节炎）
- 降低 JAK1/JAK2 抑制活性（减少副作用）

**数据准备**：
```csv
smiles,IC50_JAK1,IC50_JAK2,IC50_JAK3,Selectivity_JAK3_JAK1
# 从 BindingDB 和内部数据收集 JAK 抑制剂
# 计算选择性 = IC50_JAK1 / IC50_JAK3
```

**多靶点 SAR 分析**：
```bash
python create_sarm.py \
    --csvFile JAK_inhibitors.csv \
    --column IC50_JAK1 IC50_JAK2 IC50_JAK3 Selectivity_JAK3_JAK1 \
    --log 1 \                    # 对 IC50 取对数（更好的线性关系）
    --minimumSite1 5 \           # 提高阈值（只保留高频片段）
    --n_jobs 16 \
    --csv2excel 1
```

**SAR 树高亮选择性化合物**：
```bash
python create_tree.py \
    --fragment_core "JAK3核心骨架" \
    --workFolder JAK_SAR \
    --highlightDict "[
        {'col':'IC50_JAK3', 'type':'means', 'relation':'<', 'value':10.0},
        {'col':'Selectivity_JAK3_JAK1', 'type':'means', 'relation':'>', 'value':20.0}
    ]"
```

**关键发现**（示例）：
- **铰链区取代**：保守的氢键供体/受体（吡啶、吡唑）
- **门控区取代**：影响 JAK3 选择性的关键位点（C-5 位取代基）
- **溶剂暴露区**：可修饰以改善成药性而不影响活性

#### 2.2 TLR 抑制剂（TLR7/TLR8）

**项目包含完整案例**（`create_tree.ipynb`）：

**研究目标**：
- TLR7/TLR8 双抑制剂（治疗自身免疫疾病）
- 或 TLR7 选择性抑制剂（避免 TLR8 介导的副作用）

**数据来源**：
- 文献报道的 TLR7/8 抑制剂
- 专利化合物
- 内部合成数据

**分析流程**：
```bash
# 1. 生成 TLR7/TLR8 双靶点 SAR 表
python create_sarm.py \
    --csvFile TLR_compounds.csv \
    --column IC50_TLR7 IC50_TLR8 \
    --save_folder TLR_SAR

# 2. 构建 SAR 树（示例见 create_tree.ipynb）
# 可视化不同取代基对 TLR7/TLR8 选择性的影响
```

**关键结构特征**：
- **吲哚核心**：TLR7/8 抑制剂的常见骨架
- **吡唑并嘧啶**：增强 TLR8 选择性
- **羟甲基侧链**：影响溶解度和口服吸收

#### 2.3 ITK 抑制剂（ITK/JAK3 选择性）

**研究背景**：
- ITK（白介素-2诱导的 T 细胞激酶）
- 与 JAK3 结构相似，容易产生交叉抑制
- 目标：高 ITK 活性，低 JAK3 活性

**选择性优化策略**：
```bash
python print_sar_single_core.py \
    --csvFile ITK_series.csv \
    --core "ITK抑制剂核心骨架" \
    --actCols IC50_ITK IC50_JAK3 Selectivity_ITK_JAK3 \
    --save_folder ITK_Core_SAR
```

**分析重点**：
- 识别 ITK 选择性口袋的关键取代基
- 避免 JAK3 结合的有害基团
- 平衡活性、选择性和成药性

#### 2.4 CDK 抑制剂（CDK1/CDK2）

**完整案例研究**（`CDK1-2/` 文件夹）：

**包含内容**：
- `fragmentize_CDK2.ipynb`：CDK2 抑制剂片段化示例
- `Create_tree_AIXB7.ipynb`：CDK2 SAR 树生成
- `CDK2_SMILES/`：完整的片段化结果和 SAR 表格

**数据规模**：
- 1000+ CDK2 抑制剂（从 BindingDB 提取）
- 生成 100+ 个 2D SAR 矩阵
- 涵盖多种核心骨架（嘌呤、吡唑并嘧啶、氨基噻唑等）

**使用方法**：
```bash
# 参考 CDK1-2/fragmentize_CDK2.ipynb
# 包含完整的数据处理和分析流程
```

**应用价值**：
- **细胞周期调控**：CDK1/CDK2 在细胞分裂中的作用
- **肿瘤治疗**：CDK 抑制剂作为抗癌药物
- **选择性优化**：区分 CDK2 和其他 CDK 亚型

#### 2.5 其他激酶靶点

autoSARM 适用于所有激酶抑制剂设计：

- **EGFR/HER2**：受体酪氨酸激酶
- **BRAF/MEK**：MAPK 通路激酶
- **PI3K/AKT/mTOR**：PI3K 通路
- **SRC/ABL**：非受体酪氨酸激酶
- **Aurora A/B**：细胞有丝分裂激酶

**通用工作流程**：
1. 收集靶点抑制剂数据（BindingDB, ChEMBL, 专利）
2. 运行 `create_sarm.py` 生成 SAR 表
3. 识别保守的铰链区结合模式
4. 优化选择性口袋取代基
5. 使用 3D 对接和网格分析验证设计

---

### 3. 多靶点选择性分析

#### 3.1 双靶点选择性比值计算

**自动计算选择性**：
```python
# 在输入 CSV 中预先计算选择性
import pandas as pd
df = pd.read_csv('compounds.csv')
df['Selectivity_A_B'] = df['IC50_TargetB'] / df['IC50_TargetA']
df.to_csv('compounds_with_selectivity.csv', index=False)
```

**SAR 分析**：
```bash
python create_sarm.py \
    --csvFile compounds_with_selectivity.csv \
    --column IC50_TargetA IC50_TargetB Selectivity_A_B \
    --save_folder Selectivity_SAR
```

**高亮高选择性化合物**：
```bash
python create_tree.py \
    --highlightDict "[{'col':'Selectivity_A_B', 'type':'means', 'relation':'>', 'value':10.0}]"
```

#### 3.2 多靶点活性可视化

**三靶点示例**（如 JAK1/JAK2/JAK3）：
```csv
smiles,IC50_JAK1,IC50_JAK2,IC50_JAK3
# 在 SAR 表格中同时显示三个靶点的活性
```

**生成多目标 SAR 树**：
- 每个节点显示三个靶点的活性
- 颜色编码表示选择性模式
- 快速识别泛抑制剂 vs 选择性抑制剂

#### 3.3 选择性优化指导

**策略1：识别选择性决定基团**
```python
# 在 2D SAR 矩阵中比较：
# - 对靶点A活性高、靶点B活性低的取代基 = 靶点A选择性基团
# - 对两个靶点活性都高的取代基 = 保守结合区域
```

**策略2：结构修饰预测**
```bash
# 基于已知选择性化合物的核心
python print_sar_single_core.py \
    --core "已知选择性化合物核心" \
    --actCols IC50_TargetA IC50_TargetB Selectivity
# 分析哪些 R 基团组合可以进一步提升选择性
```

---

### 4. 专利分析与竞争情报

#### 4.1 竞争对手化合物系列分析

**数据收集**：
```python
# 从专利文献中提取化合物
# 工具：ChemDraw, MarvinSketch, OPSIN (名称转SMILES)
# 输出：patent_compounds.csv
```

**批量分析**：
```bash
python create_sarm.py \
    --csvFile patent_compounds.csv \
    --column predicted_activity \      # 可以使用预测活性
    --minimumSite1 2 \                 # 降低阈值以覆盖更多变体
    --save_folder Patent_SAR
```

**竞争分析目标**：
- 识别专利保护的核心骨架
- 列举专利权利要求中的 R 基团范围
- 发现专利未覆盖的设计空间

#### 4.2 专利权利要求覆盖范围评估

**权利要求解析**：
```
# 典型专利权利要求：
化合物通式：
  核心：A
  其中 R1 选自：甲基、乙基、丙基、...
       R2 选自：苯基、吡啶基、...
```

**覆盖范围可视化**：
```bash
python print_sar_single_core.py \
    --csvFile patent_enumerated_compounds.csv \
    --core "专利核心A" \
    --actCols predicted_IC50
# 生成 R1 × R2 矩阵，直观展示专利覆盖范围
```

#### 4.3 FTO (Freedom to Operate) 分析

**目标**：确认自己的化合物是否侵犯他人专利

**分析流程**：
```bash
# 步骤1：收集相关专利的化合物
# 步骤2：生成专利化合物的 SAR 表
python create_sarm.py --csvFile all_patent_compounds.csv

# 步骤3：检查自己的化合物是否落入专利的核心骨架和R基团范围
# 方法：子结构匹配 + SAR 表格对照
```

**决策支持**：
- ✅ 绿灯：核心骨架或 R 基团不在专利范围内
- ⚠️ 黄灯：部分重叠，需进一步法律评估
- 🛑 红灯：完全落入专利范围，需要设计绕开

---

### 5. 数据库挖掘与文献分析

#### 5.1 BindingDB、ChEMBL 等公开数据库处理

**项目包含示例**：`scripts/BindingDB_processing.ipynb`

**数据下载**：
```python
# 从 BindingDB 下载特定靶点数据
# https://www.bindingdb.org/
# 搜索：Target = "JAK3"
# 下载：TSV 格式
```

**数据清洗**：
```python
import pandas as pd

# 读取 BindingDB 数据
df = pd.read_csv('BindingDB_JAK3.tsv', sep='\t')

# 清洗步骤
df = df[df['IC50 (nM)'].notna()]           # 去除缺失值
df = df[df['Ligand SMILES'].notna()]
df['IC50_uM'] = df['IC50 (nM)'] / 1000    # 单位转换
df = df[['Ligand SMILES', 'IC50_uM']]
df.columns = ['smiles', 'IC50_uM']

# 保存清洗后的数据
df.to_csv('JAK3_clean.csv', index=False)
```

**SAR 分析**：
```bash
python create_sarm.py \
    --csvFile JAK3_clean.csv \
    --column IC50_uM \
    --save_folder JAK3_Public_SAR
```

#### 5.2 文献化合物批量提取

**工具**：
- **ChemDataExtractor**：从 PDF 提取化学结构
- **OPSIN**：化学名称转 SMILES
- **OSRA**：图片识别化学结构

**示例工作流程**：
```bash
# 步骤1：从文献 PDF 提取化合物表格
# 步骤2：使用 OPSIN 转换化学名称为 SMILES
python convert_names_to_smiles.py literature_compounds.csv

# 步骤3：SAR 分析
python create_sarm.py --csvFile literature_smiles.csv
```

#### 5.3 历史数据再分析

**场景**：项目积累了多年的合成数据，需要重新挖掘

**数据整合**：
```python
# 合并不同批次的数据
import pandas as pd
df1 = pd.read_csv('batch_2020.csv')
df2 = pd.read_csv('batch_2021.csv')
df3 = pd.read_csv('batch_2022.csv')
df_all = pd.concat([df1, df2, df3], ignore_index=True)

# 去重（相同 SMILES）
df_all = df_all.drop_duplicates(subset='smiles')
df_all.to_csv('historical_data_all.csv', index=False)
```

**重新分析**：
```bash
python create_sarm.py \
    --csvFile historical_data_all.csv \
    --column IC50_uM Selectivity Solubility LogD \
    --minimumSite1 5 \       # 更高阈值（历史数据量大）
    --n_jobs 16 \
    --csv2excel 1
```

**价值**：
- 发现之前被忽略的构效关系
- 识别历史上合成过的化合物（避免重复合成）
- 为新项目提供起点

---

## 📊 示例

### CDK2 抑制剂案例

项目包含完整的 CDK2 抑制剂分析案例（`CDK1-2/` 文件夹）：

1. **数据准备**：从 BindingDB 提取 CDK2 抑制剂数据
2. **片段化**：生成 1000+ 个片段
3. **SAR 表格**：生成 100+ 个 2D SAR 矩阵
4. **树状图**：可视化分子衍生关系

### 输出示例

**SAR 表格示例**：

```
        | R1: -CH3  | R1: -CH2CH3 | R1: -CF3
--------|-----------|-------------|----------
R2: -H  |   12.5    |    8.3      |   15.2
R2: -F  |   3.4     |    1.8      |   4.7
R2: -Cl |   5.6     |    2.9      |   6.8
```

每个单元格包含：
- 化合物结构图
- 活性统计（均值、标准差、中位数）
- 样本数量

---

## 🔬 技术特点

### 化学信息学算法

1. **分子片段化**
   - 基于环系统附近单键的智能切割
   - 环-环之间单键切割（可选）
   - 非环单键切割（可选）
   - 环系统保护
   - 芳香性保持
   - 同位素标记去除

2. **相似性计算**
   - Morgan 指纹（ECFP4，radius=2）
   - Tanimoto 相似度
   - 基于相似性的自动排序
   - 智能分子聚类

3. **子结构匹配**
   - SMARTS 模式匹配
   - 精确环系统匹配
   - 最大公共子结构（MCS）
   - 核心结构识别
   - 离子化状态容错

4. **3D 空间分析**
   - 原子坐标提取和处理
   - 3D 网格划分（可自定义单元格大小）
   - 原子密度统计（按类型和芳香性）
   - 片段空间占据评分

### 性能优化

- **并行计算**：基于 `pandarallel` 的多核并行
- **内存优化**：分块处理大数据集
- **缓存机制**：避免重复计算

### 可扩展性

- **模块化设计**：核心算法与应用分离
- **自定义插件**：支持自定义片段化规则
- **API 接口**：可集成到其他工作流

---

## 🤝 贡献指南

欢迎贡献代码、报告问题或提出建议！

### 如何贡献

1. Fork 本仓库
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 开启 Pull Request

### 报告问题

请在 [GitHub Issues](https://github.com/SongyouZhong/autoSARM/issues) 中报告问题，并提供：
- 问题描述
- 重现步骤
- 期望结果
- 实际结果
- 环境信息（Python 版本、操作系统等）

---

## 📄 许可证

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件

---

## 👥 作者

**SongyouZhong**

- GitHub: [@SongyouZhong](https://github.com/SongyouZhong)

---

## 🙏 致谢

- [RDKit](https://www.rdkit.org/) - 化学信息学核心库
- [Graphviz](https://graphviz.org/) - 图形可视化
- [Pandas](https://pandas.pydata.org/) - 数据处理

---

## 📚 相关资源

- [RDKit 文档](https://www.rdkit.org/docs/)
- [SMILES 教程](https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html)
- [药物化学基础](https://en.wikipedia.org/wiki/Medicinal_chemistry)
- [构效关系分析](https://en.wikipedia.org/wiki/Structure%E2%80%93activity_relationship)

---

## ❓ 常见问题

### Q1: 支持哪些输入格式？
A: 目前支持 CSV 格式，必须包含 SMILES 列和活性数据列。

### Q2: 可以分析多少个化合物？
A: 理论上无限制，建议单次分析 1000-10000 个化合物以获得最佳性能。

### Q3: 如何处理立体化学？
A: RDKit 会保留立体化学信息，但片段化时可能会丢失。建议预先处理立体异构体。

### Q4: 支持 3D 结构吗？
A: 支持！项目包含完整的 3D 空间分析模块（`grid_pos_utils.py`），可以：
- 处理 Schrödinger 对接结果（maegz 格式）
- 提取和分析分子的 3D 坐标
- 进行原子位置网格化分析
- 计算片段在结合口袋中的空间占据
- 生成 PyMOL 可视化

### Q5: 如何自定义片段化规则？
A: 可以修改 `utils/sarm_utils.py` 中的 `frag_mol_near_ring` 函数。通过 `pos_args` 参数控制：
- `{'RR':True, 'nRnR':False}`: 只切割环旁边的单键 + 环-环之间的单键
- `{'RR':False, 'nRnR':True}`: 切割所有非环单键
- `{'RR':True, 'nRnR':True}`: 两者结合

### Q6: 需要 Schrödinger 许可证吗？
A: 仅当需要处理对接结果时才需要 Schrödinger Suite。核心 SAR 分析功能完全独立，不依赖商业软件。

### Q7: 如何处理多个取代位点？
A: `print_sar_single_core.py` 支持自动识别核心结构上的多个虚拟原子（dummy atoms），可以处理双取代、三取代等复杂情况。

---

## 📞 联系方式

如有问题或合作意向，请通过以下方式联系：

- 📧 Email: [通过 GitHub 个人资料查看]
- 💬 GitHub Issues: [提交问题](https://github.com/SongyouZhong/autoSARM/issues)
- 🌟 Star 本项目以获取更新

---

**⭐ 如果这个项目对你有帮助，请给个 Star！**

---

*最后更新：2025年11月13日*
