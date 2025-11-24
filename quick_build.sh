#!/bin/bash

# 快速构建脚本 - 自动添加日期标签和 latest 标签

set -e

# 颜色定义
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}======================================"
echo "autoSARM 快速构建"
echo "======================================${NC}"

# 生成日期标签
DATE_TAG=$(date +%Y%m%d)
TIME_TAG=$(date +%Y%m%d-%H%M%S)

echo -e "${YELLOW}构建信息:${NC}"
echo "  日期标签: $DATE_TAG"
echo "  时间标签: $TIME_TAG"
echo ""

# 询问用户是否使用时间戳
read -p "是否使用详细时间戳 (包含时分秒)? (y/N) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    TAG=$TIME_TAG
else
    TAG=$DATE_TAG
fi

echo -e "${BLUE}开始构建...${NC}"
echo ""

# 使用 docker build 同时创建多个标签
docker build \
    -t autosarm:latest \
    -t autosarm:$TAG \
    .

echo ""
echo -e "${GREEN}✓ 构建完成！${NC}"
echo ""
echo "已创建镜像标签:"
docker images autosarm --format "  - {{.Repository}}:{{.Tag}} (大小: {{.Size}}, 创建于: {{.CreatedAt}})"
echo ""
echo -e "${BLUE}使用方法:${NC}"
echo "  # 使用 latest 标签"
echo "  docker run --rm autosarm:latest python --version"
echo ""
echo "  # 使用日期标签"
echo "  docker run --rm autosarm:$TAG python --version"
echo ""
echo "  # 使用 docker-compose"
echo "  docker-compose run --rm autosarm"
echo ""
echo "  # 启动 Jupyter"
echo "  docker-compose up autosarm-jupyter"
echo ""
