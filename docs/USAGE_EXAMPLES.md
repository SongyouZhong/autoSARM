# 📖 autoSARM 使用示例

## 🎯 快速开始

### 1. 生成 SAR 表格

#### 基本用法（使用默认参数）
```bash
python create_sarm.py \
    --csvFile data/molecules.csv \
    --column IC50_uM \
    --type smiles \
    --n_jobs 8 \
    --save_folder SAR_Results
```

#### 完整参数示例
```bash
python create_sarm.py \
    --csvFile SAR_Results/input.csv \
    --column IC50_uM pIC50 \
    --type smiles \
    --log 1 \
    --minimumSite1 3 \
    --minimumSite2 3 \
    --n_jobs 8 \
    --save_folder SAR_Results \
    --csv2excel 1
```

**参数说明**:
- `--csvFile`: **必需** - 输入CSV文件路径
- `--column`: **必需** - 活性数据列名（可以多个）
- `--type`: smiles 或 scaffold（默认：smiles）
- `--log`: 是否对数转换（0/1，默认：0）
- `--minimumSite1`: 左侧片段最小出现次数（默认：3）
- `--minimumSite2`: 右侧片段最小出现次数（默认：3）
- `--n_jobs`: 并行核心数（默认：8）
- `--save_folder`: 保存目录（默认：SAR_Tables）
- `--csv2excel`: 是否转换为Excel（0/1，默认：0）

---

### 2. 生成 SAR 树

#### 基本用法（使用默认input.csv）
```bash
python create_tree.py \
    --fragment_core "*c1ccc(*)cc1" \
    --rootTitle "Table_1_combine" \
    --workFolder ./SAR_Results \
    --maxLevel 5 \
    --treeContent "['double-cut']" \
    --highlightDict "[{'col':'IC50_uM', 'type':'mean', 'relation':'<', 'value':1.0}]"
```

#### 自定义输入文件（新功能）✨
```bash
python create_tree.py \
    --fragment_core "*c1ccc(*)cc1" \
    --rootTitle "Table_1_combine" \
    --workFolder ./SAR_Results \
    --inputFile "my_compounds.csv" \
    --maxLevel 5 \
    --treeContent "['double-cut']" \
    --highlightDict "[{'col':'IC50_uM', 'type':'mean', 'relation':'<', 'value':1.0}]"
```

**参数说明**:
- `--fragment_core`: **必需** - 核心片段SMILES（*表示连接点）
- `--rootTitle`: **必需** - 根节点表格名称
- `--workFolder`: **必需** - 工作目录（包含SAR表格）
- `--inputFile`: 输入CSV文件名（默认：`input.csv`）**新增参数** ⭐
- `--maxLevel`: 树的最大深度（默认：5）
- `--treeContent`: 树内容类型（默认：`['double-cut']`）
- `--highlightDict`: 高亮条件（JSON格式）

---

## 📂 文件结构要求

### create_sarm.py 输入文件
**必须包含列**:
- `smiles` 或 `SMILES`: 分子SMILES结构
- 活性数据列（如 `IC50_uM`）

**示例** (`input.csv`):
```csv
smiles,IC50_uM,Compound_ID
c1ccccc1,1.2,Cpd_001
CCO,2.5,Cpd_002
CCCC,0.8,Cpd_003
```

### create_tree.py 输入文件
**必须包含列**:
- `smiles` 或 `SMILES`: 分子SMILES结构
- highlightDict 中指定的活性数据列

**示例** (`my_compounds.csv`):
```csv
smiles,IC50_uM,pIC50
c1ccccc1,1.2,5.92
CCO,2.5,5.60
CCCC,0.8,6.10
```

---

## 🔍 实际应用场景

### 场景1: 处理不同来源的数据

**项目A数据**:
```bash
python create_sarm.py \
    --csvFile project_A/compounds_A.csv \
    --column IC50 \
    --save_folder SAR_Results_A

python create_tree.py \
    --workFolder ./SAR_Results_A \
    --inputFile compounds_A.csv \
    --rootTitle "Table_best" \
    --fragment_core "*c1ccccc1"
```

**项目B数据**:
```bash
python create_sarm.py \
    --csvFile project_B/compounds_B.csv \
    --column pIC50 \
    --save_folder SAR_Results_B

python create_tree.py \
    --workFolder ./SAR_Results_B \
    --inputFile compounds_B.csv \
    --rootTitle "Table_top" \
    --fragment_core "*C1CCCCC1"
```

---

### 场景2: 数据子集分析

对同一数据集的不同子集进行分析：

```bash
# 1. 生成完整SAR表格
python create_sarm.py \
    --csvFile data/all_compounds.csv \
    --column IC50_uM \
    --save_folder SAR_All

# 2. 对高活性化合物生成树
python create_tree.py \
    --workFolder ./SAR_All \
    --inputFile high_activity.csv \
    --rootTitle "Table_1" \
    --highlightDict "[{'col':'IC50_uM', 'type':'mean', 'relation':'<', 'value':0.1}]"

# 3. 对低活性化合物生成树
python create_tree.py \
    --workFolder ./SAR_All \
    --inputFile low_activity.csv \
    --rootTitle "Table_1" \
    --highlightDict "[{'col':'IC50_uM', 'type':'mean', 'relation':'>', 'value':10.0}]"
```

---

### 场景3: 多靶点分析

处理具有多个活性列的数据：

```bash
# 生成包含多个活性列的SAR表格
python create_sarm.py \
    --csvFile dual_target_data.csv \
    --column Target1_IC50 Target2_IC50 Selectivity \
    --save_folder SAR_DualTarget

# 对Target1优化
python create_tree.py \
    --workFolder ./SAR_DualTarget \
    --inputFile dual_target_data.csv \
    --rootTitle "Table_core" \
    --highlightDict "[{'col':'Target1_IC50', 'type':'mean', 'relation':'<', 'value':1.0}]"

# 对选择性优化
python create_tree.py \
    --workFolder ./SAR_DualTarget \
    --inputFile dual_target_data.csv \
    --rootTitle "Table_core" \
    --highlightDict "[{'col':'Selectivity', 'type':'mean', 'relation':'>', 'value':100}]"
```

---

## ⚙️ highlightDict 参数详解

### 基本格式
```python
"[{'col':'列名', 'type':'统计类型', 'relation':'关系', 'value':阈值}]"
```

### type 参数（支持单/复数）
- `'mean'` 或 `'means'`: 均值
- `'median'` 或 `'medians'`: 中位数
- `'std'` 或 `'stds'`: 标准差
- `'min'`: 最小值
- `'max'`: 最大值

### relation 参数
- `'<'`: 小于
- `'='`: 等于
- `'>'`: 大于

### 示例

**高亮低IC50化合物（高活性）**:
```bash
--highlightDict "[{'col':'IC50_uM', 'type':'mean', 'relation':'<', 'value':1.0}]"
```

**高亮高选择性化合物**:
```bash
--highlightDict "[{'col':'Selectivity', 'type':'median', 'relation':'>', 'value':50}]"
```

**多条件高亮**:
```bash
--highlightDict "[{'col':'IC50_uM', 'type':'mean', 'relation':'<', 'value':1.0}, {'col':'Selectivity', 'type':'median', 'relation':'>', 'value':10}]"
```

---

## 🔧 故障排查

### 问题1: FileNotFoundError - input.csv not found

**原因**: 使用默认的 `input.csv`，但文件不存在

**解决方案**: 
```bash
# 方案1: 创建 input.csv
cp your_data.csv SAR_Results/input.csv

# 方案2: 使用 --inputFile 参数
python create_tree.py \
    --inputFile "your_data.csv" \
    ...其他参数
```

---

### 问题2: KeyError: 列名不存在

**原因**: highlightDict 中的列名在输入文件中不存在

**解决方案**:
```bash
# 检查文件列名
head -1 SAR_Results/input.csv

# 确保 highlightDict 中的 col 与文件列名完全一致（区分大小写）
--highlightDict "[{'col':'正确的列名', ...}]"
```

---

### 问题3: SMILES 列不存在

**原因**: 输入文件缺少 `smiles` 或 `SMILES` 列

**解决方案**:
```bash
# 确保CSV文件包含 smiles 列（小写）
# 或在代码中添加列名映射
```

---

## 📝 最佳实践

### 1. 文件命名规范
- 使用描述性文件名: `compounds_kinase_series1.csv`
- 避免特殊字符和空格
- 使用日期标识: `compounds_2025-01-23.csv`

### 2. 参数设置建议
- **小数据集** (< 20个化合物): `--minimumSite1 2 --minimumSite2 2`
- **中等数据集** (20-100个): `--minimumSite1 3 --minimumSite2 3`（默认）
- **大数据集** (> 100个): `--minimumSite1 5 --minimumSite2 5`

### 3. 性能优化
- 根据CPU核心数调整 `--n_jobs`
- 大数据集时减小 `--maxLevel` 加快SAR树生成
- 使用 `--csv2excel 0` 跳过Excel转换以节省时间

---

## 🎯 工作流程示例

### 完整分析流程
```bash
#!/bin/bash
# 1. 准备数据
PROJECT="MyKinaseProject"
INPUT_FILE="kinase_compounds.csv"
WORK_DIR="SAR_${PROJECT}"

# 2. 生成SAR表格
python create_sarm.py \
    --csvFile data/${INPUT_FILE} \
    --column IC50_uM Selectivity \
    --type smiles \
    --log 1 \
    --n_jobs 16 \
    --save_folder ${WORK_DIR}

# 3. 查看生成的表格
echo "Generated tables:"
head -20 ${WORK_DIR}/Combine_Table_info.csv

# 4. 选择最佳表格生成树（假设 Table_10_combine 最有价值）
FRAGMENT_CORE=$(grep "Table_10_combine" ${WORK_DIR}/Combine_Table_info.csv | cut -d',' -f2)

# 5. 生成SAR树
python create_tree.py \
    --fragment_core "${FRAGMENT_CORE}" \
    --rootTitle "Table_10_combine" \
    --workFolder ./${WORK_DIR} \
    --inputFile "${INPUT_FILE}" \
    --maxLevel 4 \
    --treeContent "['double-cut']" \
    --highlightDict "[{'col':'IC50_uM', 'type':'mean', 'relation':'<', 'value':1.0}]"

# 6. 查看结果
xdg-open ${WORK_DIR}/Trees/FragTree_Table_10_combine/Table_10_combine.pdf

echo "Analysis complete! Results in ${WORK_DIR}/"
```

---

## 📚 更多资源

- **详细文档**: `TREE_GENERATION_GUIDE.md`
- **修复记录**: `FIXES_APPLIED.md`
- **错误分析**: `ERROR_ANALYSIS.md`
- **成功案例**: `SAR_TREE_SUCCESS.md`

---

**更新日期**: 2025年11月23日  
**版本**: v2.0 - 添加 --inputFile 参数支持
