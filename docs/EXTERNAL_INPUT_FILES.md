# 外部输入文件说明

本文档列出了autoSARM项目中需要用户从外部提供的所有输入文件。

## 一、必需的输入文件

### 1. CSV活性数据文件（必需）

**用途**: 包含分子SMILES和活性数据的主要输入文件

**位置**: 
- `create_sarm.py`: 通过 `--csvFile` 参数指定
- `create_tree.py`: 通过 `--inputFile` 参数指定

**代码位置**:
```python
# create_sarm.py, 第87-96行
act_csv_file = args.csvFile
dfAcat = pd.read_csv(act_csv_file)
```

**文件格式要求**:
- 格式: CSV文件
- 必须包含列:
  - `smiles`: 分子的SMILES字符串
  - 活性数据列（如 `IC50_uM`）: 通过 `--column` 参数指定

**示例**:
```csv
smiles,IC50_uM,compound_id
CCO,1.5,compound_1
c1ccccc1,2.3,compound_2
CC(C)O,0.8,compound_3
```

**命令行示例**:
```bash
python create_sarm.py --csvFile SAR_Results/input.csv --column IC50_uM --type smiles
```

---

## 二、可选的输入文件

以下文件仅在使用特定功能时需要提供。

### 2. SDF分子结构文件（可选）

**用途**: 三维分子结构文件，用于网格位置分析和片段空间定位

**位置**: 
- `utils/grid_pos_utils.py`: 通过 `sdfFolder` 参数指定文件夹

**代码位置**:
```python
# utils/grid_pos_utils.py, 第50行
sdfFiles = glob.glob(f"{sdfFolder}/*.sdf")

# utils/grid_pos_utils.py, 第28-29行
suppl = Chem.SDMolSupplier(sdf)
mol = suppl[0]
```

**文件格式要求**:
- 格式: SDF (Structure Data File)
- 来源: 分子对接软件输出（如AutoDock, Glide, GOLD等）
- 包含三维坐标信息

**使用场景**:
- 分析片段在蛋白结合口袋中的空间分布
- 计算片段位置得分
- 生成3D网格密度图

---

### 3. Maegz对接结果文件（可选）

**用途**: Schrödinger软件的分子对接结果文件

**位置**: 
- `utils/grid_pos_utils.py`: `extract_sdf()` 函数

**代码位置**:
```python
# utils/grid_pos_utils.py, 第23行
def extract_sdf(maegz, saveDir, scoreLimit):
    os.system(f"{schrodinger_root}/run {extract_sdf_fromMaegz} --maegz {maegz} ...")
```

**文件格式要求**:
- 格式: MAEGZ (Maestro compressed format)
- 来源: Schrödinger Glide或其他Schrödinger对接工具
- 需要安装Schrödinger软件

**使用场景**:
- 从Schrödinger对接结果中提取SDF文件
- 需要配合 `extract_sdf_fromMaegz.py` 脚本使用

**注意**: 需要正确配置Schrödinger路径:
```python
schrodinger_root = "/public/software/schrodinger/2021-2"
```

---

### 4. 活性数据文件（dual_target功能）（可选）

**用途**: 用于双靶点分析的活性数据

**位置**: 
- `utils/dual_target_utils.py`: `get_activity_info()` 函数

**代码位置**:
```python
# utils/dual_target_utils.py, 第154行
df_act = pd.read_csv(file_act)
df_act['matched'] = df_act.apply(lambda x: match_frag(x['Cano_SMILES'], ismarts=frag), axis=1)
```

**文件格式要求**:
- 格式: CSV文件
- 必须包含列:
  - `Cano_SMILES`: 规范化的SMILES
  - 活性数据列

**使用场景**:
- 双靶点SAR分析
- 片段-活性关联分析

---

### 5. SAR表格信息文件（tree构建）（可选）

**用途**: 用于构建SAR树的输入数据

**位置**: 
- `create_tree.py`: 从 `workFolder` 读取

**代码位置**:
```python
# create_tree.py, 第62行
df_tmp = pd.read_csv(f"{workFolder}/Combine_Table_info.csv")

# create_tree.py, 第69行
df_tmp = pd.read_csv(f"{workFolder}/singleCut_Table_info.csv")
```

**文件格式要求**:
- 格式: CSV文件
- 这些文件通常由 `create_sarm.py` 生成
- 必须包含列:
  - `Key2` 或 `SMILES`: 片段的SMILES

**使用场景**:
- 基于之前生成的SAR表构建片段树
- 层次化分析SAR关系

**注意**: 这些文件虽然是外部输入（对于`create_tree.py`），但通常是由`create_sarm.py`程序生成的。

---

## 三、文件准备指南

### 准备CSV活性数据文件

1. **收集数据**: 从实验数据、文献或数据库中收集化合物SMILES和活性数据
2. **格式化**: 确保CSV文件包含必需的列
3. **验证SMILES**: 使用RDKit或其他工具验证SMILES的有效性
4. **保存**: 保存为UTF-8编码的CSV文件

### 准备SDF文件（如需要）

1. **分子对接**: 使用对接软件（AutoDock Vina, Glide等）进行分子对接
2. **导出结果**: 将对接姿态导出为SDF格式
3. **组织文件**: 将所有SDF文件放在同一文件夹中
4. **命名规范**: 建议使用化合物ID作为文件名

### 从Maegz提取SDF（如需要）

```bash
# 需要Schrödinger软件和extract_sdf_fromMaegz.py脚本
python extract_sdf.py --maegz docking_results.maegz --saveDir ./sdf_output --scoreLimit 6.0
```

---

## 四、常见问题

### Q1: CSV文件路径不存在怎么办？
**A**: 确保使用绝对路径或相对于工作目录的正确路径。可以添加文件存在性检查：
```python
import os
if not os.path.exists(csv_file):
    raise FileNotFoundError(f"CSV file not found: {csv_file}")
```

### Q2: SMILES格式错误怎么办？
**A**: 程序会自动过滤无效的SMILES。查看日志文件了解哪些分子被过滤。

### Q3: SDF文件太大怎么办？
**A**: 
- 只提取对接得分较好的姿态
- 使用 `scoreLimit` 参数过滤
- 分批处理大型数据集

### Q4: 不同格式的分子文件如何转换？
**A**: 可以使用以下工具：
- Open Babel: `obabel input.mol2 -O output.sdf`
- RDKit: 程序化转换
- Schrödinger: 使用 `structconvert` 命令

---

## 五、文件示例

### 示例1: 基本CSV输入文件

```csv
smiles,IC50_uM,compound_id,molecular_weight
CCO,1.5,CHEMBL123,46.07
c1ccccc1,2.3,CHEMBL456,78.11
CC(C)O,0.8,CHEMBL789,60.10
CCN(CC)CC,3.2,CHEMBL321,101.19
```

### 示例2: 多列活性数据

```csv
smiles,IC50_Target1,IC50_Target2,Solubility,LogP
CCO,1.5,3.2,high,-0.31
c1ccccc1,2.3,1.8,medium,2.13
CC(C)O,0.8,2.1,high,0.05
```

### 示例3: 带有更多属性的输入文件

```csv
smiles,pIC50,Ki_nM,compound_id,batch,date,notes
CCO,8.82,15.2,COMP001,Batch_A,2024-01-15,Active
c1ccccc1,8.64,23.1,COMP002,Batch_A,2024-01-15,Reference
CC(C)O,9.10,8.0,COMP003,Batch_B,2024-01-20,Lead
```

---

## 六、总结

### 必需文件（1个）
1. ✅ **CSV活性数据文件** - 包含SMILES和活性数据

### 可选文件（根据功能需求）
2. ⭕ **SDF文件** - 用于3D空间分析
3. ⭕ **Maegz文件** - Schrödinger对接结果
4. ⭕ **活性数据文件** - 双靶点分析
5. ⭕ **SAR表格文件** - 树构建功能

### 建议工作流程

```
1. 准备CSV活性数据
   ↓
2. 运行 create_sarm.py 生成SAR表
   ↓
3. （可选）运行 create_tree.py 构建SAR树
   ↓
4. （可选）如需3D分析，准备SDF文件
   ↓
5. （可选）使用grid_pos功能进行空间分析
```

---

**最后更新**: 2025年11月23日
**维护者**: autoSARM项目组
