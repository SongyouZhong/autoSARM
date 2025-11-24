#!/bin/bash

# autoSARM Docker 构建和测试脚本

set -e  # 遇到错误立即退出

echo "======================================"
echo "autoSARM Docker 构建脚本"
echo "======================================"

# 颜色定义
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 函数：打印成功消息
print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

# 函数：打印错误消息
print_error() {
    echo -e "${RED}✗ $1${NC}"
}

# 函数：打印警告消息
print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

# 检查 Docker 是否安装
check_docker() {
    echo "检查 Docker 安装..."
    if ! command -v docker &> /dev/null; then
        print_error "Docker 未安装，请先安装 Docker"
        exit 1
    fi
    print_success "Docker 已安装: $(docker --version)"
}

# 检查 Docker Compose 是否安装
check_docker_compose() {
    echo "检查 Docker Compose 安装..."
    if ! command -v docker-compose &> /dev/null; then
        print_warning "Docker Compose 未安装，将使用原生 Docker 命令"
        return 1
    fi
    print_success "Docker Compose 已安装: $(docker-compose --version)"
    return 0
}

# 构建 Docker 镜像
build_image() {
    echo ""
    echo "======================================"
    echo "开始构建 Docker 镜像..."
    echo "======================================"
    
    # 生成日期标签 (格式: YYYYMMDD)
    DATE_TAG=$(date +%Y%m%d)
    echo "日期标签: $DATE_TAG"
    
    if check_docker_compose; then
        # 使用 docker-compose 构建
        docker-compose build
        
        # 为镜像添加日期标签
        echo "添加日期标签: autosarm:$DATE_TAG"
        docker tag autosarm:latest autosarm:$DATE_TAG
    else
        # 直接使用 docker build，同时创建 latest 和日期标签
        docker build -t autosarm:latest -t autosarm:$DATE_TAG .
    fi
    
    print_success "镜像构建完成"
    echo "可用标签:"
    echo "  - autosarm:latest"
    echo "  - autosarm:$DATE_TAG"
}

# 测试镜像
test_image() {
    echo ""
    echo "======================================"
    echo "测试 Docker 镜像..."
    echo "======================================"
    
    echo "1. 测试 Python 环境..."
    docker run --rm autosarm:latest python --version
    print_success "Python 环境正常"
    
    echo ""
    echo "2. 测试 RDKit 安装..."
    docker run --rm autosarm:latest python -c "from rdkit import Chem; print('RDKit version:', Chem.rdMolDescriptors.Properties.GetRDKitVersion())"
    print_success "RDKit 安装正常"
    
    echo ""
    echo "3. 测试 Pandas 安装..."
    docker run --rm autosarm:latest python -c "import pandas as pd; print('Pandas version:', pd.__version__)"
    print_success "Pandas 安装正常"
    
    echo ""
    echo "4. 测试 Graphviz 安装..."
    docker run --rm autosarm:latest python -c "import graphviz; print('Graphviz 可用')"
    print_success "Graphviz 安装正常"
    
    echo ""
    echo "5. 列出已安装的包..."
    docker run --rm autosarm:latest micromamba list
}

# 运行示例
run_example() {
    echo ""
    echo "======================================"
    echo "运行示例（如果有测试数据）..."
    echo "======================================"
    
    if [ -f "SAR_Results/input.csv" ]; then
        print_warning "运行 SAR 分析示例..."
        docker run --rm \
            -v $(pwd):/app \
            autosarm:latest \
            python create_sarm.py \
                --csvFile SAR_Results/input.csv \
                --column IC50_uM \
                --type smiles \
                --log 1 \
                --minimumSite1 3 \
                --minimumSite2 3 \
                --n_jobs 4 \
                --save_folder SAR_Results_test
        print_success "示例运行成功"
    else
        print_warning "未找到测试数据文件 SAR_Results/input.csv，跳过示例运行"
    fi
}

# 清理旧镜像和容器
cleanup() {
    echo ""
    echo "======================================"
    echo "清理旧的 Docker 镜像和容器..."
    echo "======================================"
    
    # 停止运行的容器
    docker stop autosarm 2>/dev/null || true
    docker stop autosarm-jupyter 2>/dev/null || true
    
    # 删除旧容器
    docker rm autosarm 2>/dev/null || true
    docker rm autosarm-jupyter 2>/dev/null || true
    
    # 清理未使用的镜像（可选）
    read -p "是否清理未使用的 Docker 镜像？(y/N) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        docker image prune -f
        print_success "清理完成"
    fi
}

# 显示使用帮助
show_help() {
    cat << EOF
使用方法: ./build.sh [选项]

选项:
    build       构建 Docker 镜像
    test        测试 Docker 镜像
    example     运行示例
    cleanup     清理旧的镜像和容器
    all         执行所有步骤（构建、测试、运行示例）
    help        显示此帮助信息

示例:
    ./build.sh build        # 只构建镜像
    ./build.sh test         # 只测试镜像
    ./build.sh all          # 构建并测试
EOF
}

# 主函数
main() {
    check_docker
    
    case "${1:-all}" in
        build)
            build_image
            ;;
        test)
            test_image
            ;;
        example)
            run_example
            ;;
        cleanup)
            cleanup
            ;;
        all)
            build_image
            test_image
            run_example
            ;;
        help|--help|-h)
            show_help
            ;;
        *)
            print_error "未知选项: $1"
            show_help
            exit 1
            ;;
    esac
    
    echo ""
    echo "======================================"
    print_success "完成！"
    echo "======================================"
    echo ""
    echo "可用镜像:"
    docker images autosarm --format "  - {{.Repository}}:{{.Tag}} (创建于: {{.CreatedAt}})"
    echo ""
    echo "接下来可以："
    echo "  1. 运行交互式 shell:    docker-compose run --rm autosarm"
    echo "  2. 启动 Jupyter:        docker-compose up autosarm-jupyter"
    echo "  3. 查看使用说明:        cat DOCKER_USAGE.md"
    echo "  4. 查看所有标签:        docker images autosarm"
}

# 执行主函数
main "$@"
