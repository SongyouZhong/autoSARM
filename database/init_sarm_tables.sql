-- ============================================================
-- AutoSARM 任务参数表
-- 依赖 AstraMolecula 主数据库的 tasks 表
-- ============================================================

-- sarm_task_params - SARM 分析任务参数表
-- task_subtype 区分 'sarm' (矩阵生成) 和 'tree' (SAR 树生成)
CREATE TABLE IF NOT EXISTS sarm_task_params (
  id CHAR(32) NOT NULL,
  task_id CHAR(36) NOT NULL,
  task_subtype VARCHAR(20) NOT NULL DEFAULT 'sarm',

  -- SARM 矩阵生成参数 (task_subtype = 'sarm')
  csv_filename VARCHAR(255) DEFAULT 'compounds.csv',
  analysis_type VARCHAR(20) DEFAULT 'smiles',          -- 'smiles' 或 'scaffold'
  value_columns TEXT DEFAULT '[]',                      -- JSON 数组，如 '["IC50", "Ki"]'
  log_transform BOOLEAN DEFAULT FALSE,
  minimum_site1 DECIMAL(10,2) DEFAULT 3,
  minimum_site2 DECIMAL(10,2) DEFAULT 3,
  n_jobs INT DEFAULT 8,
  csv2excel BOOLEAN DEFAULT FALSE,

  -- SAR 树生成参数 (task_subtype = 'tree')
  fragment_core VARCHAR(1024) DEFAULT NULL,
  root_title VARCHAR(255) DEFAULT NULL,
  input_file VARCHAR(255) DEFAULT 'input.csv',
  tree_content TEXT DEFAULT '["double-cut"]',           -- JSON 数组
  highlight_dict TEXT DEFAULT '[]',                     -- JSON 数组
  max_level INT DEFAULT 5,

  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

  PRIMARY KEY (id),
  FOREIGN KEY (task_id) REFERENCES tasks(id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_sarm_task_params_task_id ON sarm_task_params(task_id);
CREATE INDEX IF NOT EXISTS idx_sarm_task_params_subtype ON sarm_task_params(task_subtype);

-- ============================================================
-- Done
-- ============================================================
