# Docker 使用说明

本文档介绍如何使用 Docker 和 micromamba 运行 autoSARM 项目。

## 📋 前置要求

- Docker Engine 20.10+ 
- Docker Compose 1.29+ (可选，推荐)
- 至少 4GB 可用磁盘空间

## 🚀 快速开始

### 方法 1：使用 Docker Compose（推荐）

#### 1.1 构建并启动容器

```bash

#快速构建（推荐）
./quick_build.sh

# 构建镜像
docker-compose build

# 启动交互式 Shell
docker-compose run --rm autosarm

# 或者启动 Jupyter Notebook
docker-compose up autosarm-jupyter
```

#### 1.2 在容器中运行分析

进入容器后，可以运行以下命令：

```bash
# 生成 SAR 表格
python create_sarm.py \
    --csvFile SAR_Results/input.csv \
    --column IC50_uM \
    --type smiles \
    --log 1 \
    --minimumSite1 3 \
    --minimumSite2 3 \
    --n_jobs 8 \
    --save_folder SAR_Results \
    --csv2excel 1

# 生成 SAR 树
python create_tree.py \
    --fragment_core "*CN1CCC(c2ccc3[nH]c(-c4cc(CO*)c5ncnn5c4)c(C(C)C)c3c2)CC1" \
    --rootTitle "Table_100_combine" \
    --workFolder ./SAR_Results \
    --maxLevel 5
```

### 方法 2：使用原生 Docker 命令

#### 2.1 构建镜像

```bash
docker build -t autosarm:latest .
```

#### 2.2 运行容器

```bash
# 交互式模式
docker run -it --rm \
    -v $(pwd):/app \
    -v $(pwd)/SAR_Results:/app/SAR_Results \
    autosarm:latest

# 运行特定命令
docker run --rm \
    -v $(pwd):/app \
    autosarm:latest \
    python create_sarm.py --csvFile SAR_Results/input.csv --column IC50_uM --type smiles
```

#### 2.3 运行 Jupyter Notebook

```bash
docker run -it --rm \
    -p 8888:8888 \
    -v $(pwd):/app \
    autosarm:latest \
    jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser --allow-root
```

然后在浏览器中打开 `http://localhost:8888`

## 📁 数据卷挂载说明

容器中的重要目录映射：

| 容器内路径 | 宿主机路径 | 说明 |
|-----------|-----------|------|
| `/app` | `.` (项目根目录) | 项目代码和配置文件 |
| `/app/data` | `./data` | 输入数据文件 |
| `/app/SAR_Results` | `./SAR_Results` | SAR 分析结果 |

## 🔧 环境配置

### 查看已安装的包

```bash
# 进入容器
docker-compose run --rm autosarm bash

# 查看 conda 环境
micromamba list

# 查看 Python 版本
python --version
```

### 安装额外的包

如果需要安装额外的 Python 包：

```bash
# 使用 micromamba
micromamba install -n base -c conda-forge package_name

# 或使用 pip
pip install package_name
```

**注意**：容器重启后，未保存到镜像的包会丢失。如需永久安装，需修改 `env.yaml` 并重新构建镜像。

## 🐛 常见问题

### 1. 权限问题

如果遇到文件权限问题，可以在运行容器时添加用户映射：

```bash
docker run -it --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/app \
    autosarm:latest
```

### 2. 内存不足

如果处理大型数据集时内存不足，可以增加 Docker 内存限制：

```bash
docker run -it --rm \
    --memory="8g" \
    -v $(pwd):/app \
    autosarm:latest
```

### 3. Graphviz 相关错误

如果遇到 Graphviz 错误，确保系统级 Graphviz 已安装。可以重新构建镜像：

```bash
docker-compose build --no-cache
```

### 4. RDKit 导入错误

确保使用的是容器内的 Python 环境：

```bash
# 在容器内检查
which python
# 应该输出: /opt/conda/bin/python
```

## 📝 开发建议

### 1. 使用开发模式

在 `docker-compose.yml` 中已经配置了代码目录挂载，修改代码后无需重新构建镜像。

### 2. 调试 Python 代码

```bash
# 启动容器并进入 IPython
docker-compose run --rm autosarm ipython

# 或使用 Jupyter
docker-compose up autosarm-jupyter
```

### 3. 运行测试

```bash
docker-compose run --rm autosarm pytest tests/
```

## 🔒 生产环境部署

### 1. 多阶段构建优化（可选）

可以创建一个优化的生产版 Dockerfile：

```dockerfile
# 见 Dockerfile.prod
```

### 2. 安全性建议

- 不要在生产环境中使用 `--allow-root` 运行 Jupyter
- 设置 Jupyter token 或密码
- 使用非 root 用户运行容器
- 定期更新基础镜像

## 📚 更多信息

- [Docker 官方文档](https://docs.docker.com/)
- [Micromamba 文档](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
- [autoSARM 项目 README](./README.md)

## 🤝 支持

如果遇到问题，请查看：
1. 项目 README.md
2. GitHub Issues
3. Docker logs: `docker-compose logs`
