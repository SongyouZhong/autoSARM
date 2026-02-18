# AutoSARM 用户参数说明

本文档说明通过 AstraMolecula 前端（或 API）提交 AutoSARM 任务时，**用户需要/可以配置的参数**。

---

## 1. 通用参数（所有任务类型共用）

这些字段写入主表 `tasks`，由前端/API 在创建任务时生成：

| 参数 | 类型 | 是否必填 | 说明 |
|------|------|----------|------|
| `task_type` | string | **必填** | 固定值 `sarm_analysis` |
| `task_subtype` | string | **必填** | `sarm`（SARM 矩阵生成）或 `tree`（SAR 树生成） |

> `task_type` 写入 `tasks` 表；`task_subtype` 写入 `sarm_task_params` 表。

---

## 2. SARM 矩阵生成参数（`task_subtype = 'sarm'`）

### 2.1 用户必须提供

| 参数 | 数据库字段 | 类型 | 说明 |
|------|-----------|------|------|
| **化合物 CSV 文件** | — | 文件上传 | 包含 SMILES 和活性数据的 CSV 文件，上传到 SeaweedFS 的 `jobs/sarm_analysis/{job_id}/input/` 路径下 |
| **活性列名称** | `value_columns` | JSON 数组 | CSV 中用于 SAR 分析的数值活性列名。例如 `["IC50_uM"]` 或 `["IC50", "Ki"]`。支持多列 |

### 2.2 用户可选配置（有默认值）

| 参数 | 数据库字段 | 类型 | 默认值 | 说明 |
|------|-----------|------|--------|------|
| CSV 文件名 | `csv_filename` | string | `compounds.csv` | 上传的 CSV 文件名（若前端固定命名则用户无需关心） |
| 分析类型 | `analysis_type` | string | `smiles` | `smiles` — 基于完整 SMILES 分析；`scaffold` — 基于 Murcko 骨架分析 |
| 对数变换 | `log_transform` | boolean | `false` | 是否对活性值取对数（ln）。当活性值跨越多个数量级时建议开启 |
| Site1 最小计数 | `minimum_site1` | number | `3` | 第一切割位点片段的最小出现次数。低于此值的片段将被过滤 |
| Site2 最小计数 | `minimum_site2` | number | `3` | 第二切割位点片段的最小出现次数。低于此值的片段将被过滤 |
| 并行任务数 | `n_jobs` | integer | 自动检测 | 分子片段化的并行进程数。默认使用容器可用 CPU 的 80%，一般无需用户设置 |
| 导出 Excel | `csv2excel` | boolean | `false` | 是否将结果同时导出为带分子结构图的 Excel 文件。开启后处理时间会明显增加 |

### 2.3 CSV 文件格式要求

用户上传的 CSV 文件需满足以下条件：

```
smiles,IC50_uM,Ki_nM
CC(=O)Oc1ccccc1C(=O)O,3.2,15.6
CC(C)Cc1ccc(cc1)C(C)C(=O)O,0.8,4.2
...
```

- **必须包含 `smiles` 列**（不区分大小写），内容为有效的 SMILES 字符串
- **必须包含至少一个数值活性列**，列名与 `value_columns` 参数对应
- 支持标准 CSV 格式（逗号分隔）
- 建议去除重复和无效 SMILES 行（系统也会自动处理）

---

## 3. SAR 树生成参数（`task_subtype = 'tree'`）

SAR 树通常在 SARM 矩阵生成完成后使用，依赖先前的分析结果。

### 3.1 用户必须提供

| 参数 | 数据库字段 | 类型 | 说明 |
|------|-----------|------|------|
| **片段核心** | `fragment_core` | string | 用于构建衍生树的核心片段 SMARTS/SMILES，如 `c1ccc(cc1)C(=O)` |
| **根节点标题** | `root_title` | string | SAR 树根节点的显示名称，如 `Benzamide core` |

### 3.2 用户可选配置（有默认值）

| 参数 | 数据库字段 | 类型 | 默认值 | 说明 |
|------|-----------|------|--------|------|
| 输入文件 | `input_file` | string | `input.csv` | SAR 树输入数据文件名 |
| 树内容类型 | `tree_content` | JSON 数组 | `["double-cut"]` | 树中展示的内容类型。可选值：`double-cut`、`single-cut` 或组合 |
| 高亮配置 | `highlight_dict` | JSON 数组 | `[]` | 需要在树中高亮显示的化合物/片段列表 |
| 最大层级 | `max_level` | integer | `5` | 树的最大展开层数。层数越多计算越慢，建议 3-7 |

---

## 4. 文件存储路径约定

前端上传文件时需遵循以下 SeaweedFS 路径规范：

```
jobs/sarm_analysis/{job_id}/
├── input/
│   ├── compounds.csv          ← 用户上传的 CSV 文件
│   └── input.csv              ← Tree 任务的输入文件（可选）
└── output/                    ← Worker 自动生成，用户无需操作
    └── SAR_Results/
        ├── Left_Table/
        ├── Right_Table/
        ├── Combine_Table/
        ├── Left_Table_info.csv
        ├── Right_Table_info.csv
        ├── Combine_Table_info.csv
        └── singleCut_Table_info.csv
```

---

## 5. 任务创建流程（前端/API 视角）

创建一个 SARM 矩阵分析任务需要执行三步操作：

### Step 1：上传文件到 SeaweedFS

```
POST {filer_endpoint}/buckets/astramolecula/jobs/sarm_analysis/{job_id}/input/compounds.csv
Content-Type: multipart/form-data
```

### Step 2：在 `tasks` 表中创建任务记录

```sql
INSERT INTO tasks (id, user_id, task_type, job_dir, status)
VALUES (
  '{job_id}',                           -- UUID
  '{user_id}',                          -- 用户ID
  'sarm_analysis',                      -- 固定值
  'jobs/sarm_analysis/{job_id}',        -- SeaweedFS 路径前缀
  'pending'                             -- Worker 会自动获取
);
```

### Step 3：在 `sarm_task_params` 表中写入任务参数

```sql
INSERT INTO sarm_task_params (id, task_id, task_subtype, csv_filename, value_columns, analysis_type, log_transform, minimum_site1, minimum_site2, csv2excel)
VALUES (
  '{param_id}',                         -- 32 位唯一ID
  '{job_id}',                           -- 关联的任务ID
  'sarm',                               -- 子类型
  'compounds.csv',                      -- CSV 文件名
  '["IC50_uM"]',                        -- 活性列（JSON 数组）
  'smiles',                             -- 分析类型
  false,                                -- 是否对数变换
  3,                                    -- Site1 最小计数
  3,                                    -- Site2 最小计数
  false                                 -- 是否导出 Excel
);
```

任务创建后，空闲的 Worker 实例将自动轮询获取 `pending` 状态的任务并开始处理。

---

## 6. 任务状态查询

用户可通过查询 `tasks` 表获取任务状态：

| 状态 | 含义 |
|------|------|
| `pending` | 等待中，尚未被 Worker 领取 |
| `processing` | 正在处理中 |
| `finished` | 处理完成，结果已上传到 SeaweedFS |
| `failed` | 处理失败（可查看 Worker 日志定位原因） |
| `cancelled` | 已取消 |

结果文件在任务 `finished` 后可从 SeaweedFS 下载：

```
GET {filer_endpoint}/buckets/astramolecula/jobs/sarm_analysis/{job_id}/output/SAR_Results/
```

---

## 7. 运维参数（非用户参数）

以下参数由运维人员在部署时通过 `.env` 文件或 Docker 环境变量配置，**用户无需关心**：

| 参数 | 环境变量 | 默认值 | 说明 |
|------|---------|--------|------|
| 数据库地址 | `DB_HOST` | `localhost` | PostgreSQL 主机 |
| 数据库端口 | `DB_PORT` | `5432` | PostgreSQL 端口 |
| 数据库用户 | `DB_USER` | `admin` | 数据库用户名 |
| 数据库密码 | `DB_PASSWORD` | `secret` | 数据库密码 |
| 数据库名 | `DB_NAME` | `mydatabase` | 数据库名称 |
| SeaweedFS 地址 | `SEAWEED_FILER_ENDPOINT` | `http://localhost:8888` | Filer API 地址 |
| 存储桶 | `SEAWEED_BUCKET` | `astramolecula` | SeaweedFS bucket |
| 轮询间隔 | `POLL_INTERVAL` | `30` | Worker 检查新任务的间隔（秒） |
| 日志级别 | `LOG_LEVEL` | `INFO` | 日志详细程度 |
| 临时目录 | `TEMP_DIR` | `/tmp/autosarm` | 计算临时文件目录 |
| Worker 端口 | `SARM_PORT_1/2/3` | `8030/8031/8032` | 各 Worker 健康检查端口映射 |
