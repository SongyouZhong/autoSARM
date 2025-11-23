# 🌲 SAR树生成成功报告

**生成时间**: 2025年11月23日 21:45  
**状态**: ✅ 完成

---

## 📊 生成结果

### 输出文件

```
SAR_Results/Trees/FragTree_Table_1_combine/
├── Table_1_combine.pdf (1.3 MB)          # PDF格式的SAR树可视化
├── Table_1_combine                        # GraphViz源文件
├── Combine_Table_info_Tree_Table_1_combine.txt (1.1 MB)  # 文本格式的树结构
└── Images/ (143张图片)                     # 分子结构图像
```

### 统计数据

| 项目 | 数值 |
|------|------|
| **根节点** | `*c1ccc(*)cc1` (苯环双取代) |
| **最大深度** | 3层 |
| **生成的分子图像** | **143张** |
| **PDF文件大小** | 1.3 MB |
| **文本树大小** | 1.1 MB |

---

## 🔧 修复的问题

### 1. 文件路径问题
**问题**: `create_tree.py` 期望文件在 `workFolder/SAR_Tables/` 下，但实际在 `workFolder/` 下

**修复**:
```python
# 修改前
df_tmp=pd.read_csv(f"{workFolder}/SAR_Tables/Combine_Table_info.csv")

# 修改后
df_tmp=pd.read_csv(f"{workFolder}/Combine_Table_info.csv")
```

**影响文件**: `create_tree.py` (2处修改)

### 2. highlightDict type参数兼容性
**问题**: `KeyError: 'means'` - 代码期望 `'mean'` 但用户使用 `'means'`

**修复**: 添加复数形式支持
```python
# 支持 means/mean, medians/median, stds/std
if itype == 'means':
    itype = 'mean'
elif itype == 'medians':
    itype = 'median'
elif itype == 'stds':
    itype = 'std'
```

**影响文件**: `utils/tree_utils.py`

### 3. MCS模块弃用（已在之前修复）
```python
# utils/grid_pos_utils.py
from rdkit.Chem import rdFMCS as MCS  # ✅ 已修复
```

---

## 📖 使用说明

### 查看生成的SAR树

**方法1**: PDF查看器
```bash
# 使用系统默认PDF查看器
xdg-open SAR_Results/Trees/FragTree_Table_1_combine/Table_1_combine.pdf

# 或使用其他PDF工具
evince SAR_Results/Trees/FragTree_Table_1_combine/Table_1_combine.pdf
okular SAR_Results/Trees/FragTree_Table_1_combine/Table_1_combine.pdf
```

**方法2**: 查看文本树结构
```bash
cat SAR_Results/Trees/FragTree_Table_1_combine/Combine_Table_info_Tree_Table_1_combine.txt | less
```

**方法3**: 浏览分子图像
```bash
ls SAR_Results/Trees/FragTree_Table_1_combine/Images/
# 查看单个图像
display SAR_Results/Trees/FragTree_Table_1_combine/Images/L0_0.png
```

---

## 🎨 树结构说明

### 层级结构

```
Level 0 (根节点):
  *c1ccc(*)cc1  (苯环双取代位点)
  
Level 1 (主要分支):
  ├─ *Nc1ccc(*)cc1  (苯胺取代, 26个化合物)
  ├─ *Oc1ccc(*)cc1  (苯酚取代, 18个化合物)
  ├─ *c1cc2c(*)ncnc2cc1F  (氟取代喹唑啉, 12个化合物)
  ├─ *c1cc2c(*)ncnc2cc1O  (羟基喹唑啉, 10个化合物)
  └─ *c1cc2ncnc(Nc3ccccc3)c2cc1*  (苯胺喹唑啉, 9个化合物)

Level 2-3:
  更细分的取代基组合
```

### 高亮规则

根据 `highlightDict` 设置：
- **列**: IC50_uM
- **类型**: 均值 (means)
- **条件**: < 1.0 μM
- **含义**: 高活性化合物会被突出显示

---

## 🔍 树分析建议

### 1. 识别高活性片段
查找树中高亮的节点，这些代表平均IC50 < 1.0 μM的片段组合

### 2. 活性差异比较
比较同一层级不同分支的活性差异，识别有利/不利取代基

### 3. 优化路径
从根节点到高活性叶节点的路径代表了优化方向

### 4. 数据量评估
节点中的 `(m, n)count` 表示SAR矩阵大小和化合物数量
- 数量越多，SAR趋势越可靠

---

## ⚙️ 生成参数

```bash
python create_tree.py \
    --fragment_core "*c1ccc(*)cc1" \
    --rootTitle "Table_1_combine" \
    --workFolder ./SAR_Results \
    --maxLevel 3 \
    --treeContent "['double-cut']" \
    --highlightDict "[{'col':'IC50_uM', 'type':'means', 'relation':'<', 'value':1.0}]"
```

### 参数说明

| 参数 | 值 | 说明 |
|------|-----|------|
| `fragment_core` | `*c1ccc(*)cc1` | 苯环双取代核心 |
| `rootTitle` | `Table_1_combine` | 根表格名称 |
| `workFolder` | `./SAR_Results` | 数据目录 |
| `maxLevel` | `3` | 树的最大深度 |
| `treeContent` | `['double-cut']` | 使用双切割表格 |
| `highlightDict.col` | `IC50_uM` | 高亮依据的列 |
| `highlightDict.type` | `means` | 使用均值 |
| `highlightDict.relation` | `<` | 小于阈值 |
| `highlightDict.value` | `1.0` | 阈值1.0 μM |

---

## 📝 生成其他表格的SAR树

### 步骤1: 查看可用表格
```bash
head -20 SAR_Results/Combine_Table_info.csv
```

### 步骤2: 选择感兴趣的表格
找到 `Items_count` 较大的表格（数据点多，SAR更可靠）

### 步骤3: 获取fragment_core
从 `Combine_Table_info.csv` 的 `Key2` 列获取片段SMILES

### 步骤4: 生成树
```bash
python create_tree.py \
    --fragment_core "YOUR_FRAGMENT" \
    --rootTitle "Table_XXX_combine" \
    --workFolder ./SAR_Results \
    --maxLevel 5 \
    --treeContent "['double-cut']" \
    --highlightDict "[{'col':'IC50_uM', 'type':'mean', 'relation':'<', 'value':1.0}]"
```

**注意**: `type` 参数现在支持：
- 单数形式: `mean`, `median`, `std`
- 复数形式: `means`, `medians`, `stds` (都可以)

---

## ⚠️ 已知问题

### 字体大小警告
```
The new font size 640 is above the current maximum (40).
```
**影响**: 仅警告，不影响功能  
**原因**: RDKit绘图字体大小设置  
**解决**: 可以忽略，或在 `tree_utils.py` 中调整字体大小参数

---

## ✅ 完成清单

- ✅ 修复 `create_tree.py` 文件路径问题
- ✅ 修复 `tree_utils.py` highlightDict type参数兼容性
- ✅ 修复 `grid_pos_utils.py` MCS模块弃用
- ✅ 成功生成SAR树 (143个节点图像)
- ✅ 生成PDF可视化文件 (1.3 MB)
- ✅ 生成文本树结构 (1.1 MB)

---

## 🎯 下一步建议

1. **分析PDF树**: 打开 `Table_1_combine.pdf` 查看完整的SAR树
2. **识别活性模式**: 找出高亮节点的共同特征
3. **生成更多树**: 对其他有价值的表格生成SAR树
4. **导出数据**: 从树中提取高活性片段用于虚拟筛选

---

**🎉 SAR树生成完成！所有代码问题已修复，系统运行正常。**

*工具: autoSARM*  
*环境: autoSAR (Python 3.10 + RDKit 2025.03.6)*  
*修复时间: 2025年11月23日*
