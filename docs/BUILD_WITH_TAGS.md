# Docker 镜像标签构建说明

## 概述

项目现已支持自动添加日期标签，同时保留 `latest` 标签。这样可以：
- 使用 `latest` 标签始终获取最新构建
- 使用日期标签追踪和回退到特定版本
- 方便版本管理和问题排查

## 镜像标签格式

构建时会自动创建以下标签：

1. **latest** - 始终指向最新构建
   ```bash
   autosarm:latest
   ```

2. **日期标签** - 格式：YYYYMMDD
   ```bash
   autosarm:20250124
   ```

3. **时间戳标签** (可选) - 格式：YYYYMMDD-HHMMSS
   ```bash
   autosarm:20250124-143052
   ```

## 构建方法

### 方法 1: 使用快速构建脚本（推荐）

```bash
./quick_build.sh
```

这个脚本会：
- 自动生成日期标签
- 询问是否使用详细时间戳
- 同时创建 `latest` 和日期标签
- 显示所有可用镜像

### 方法 2: 使用完整构建脚本

```bash
./build.sh build
```

这个脚本会：
- 使用 docker-compose 构建
- 自动添加日期标签
- 运行测试验证

### 方法 3: 手动构建

```bash
# 获取当前日期
DATE_TAG=$(date +%Y%m%d)

# 构建并同时创建多个标签
docker build -t autosarm:latest -t autosarm:$DATE_TAG .
```

### 方法 4: 使用 docker-compose

```bash
# 设置环境变量（可选）
export IMAGE_TAG=$(date +%Y%m%d)

# 构建
docker-compose build

# 添加日期标签
docker tag autosarm:latest autosarm:$IMAGE_TAG
```

## 使用镜像

### 使用 latest 标签

```bash
# 运行命令
docker run --rm autosarm:latest python --version

# 交互式 shell
docker run -it --rm autosarm:latest /bin/bash

# 使用 docker-compose
docker-compose run --rm autosarm
```

### 使用特定日期标签

```bash
# 运行特定日期的镜像
docker run --rm autosarm:20250124 python create_sarm.py --help

# 回退到特定版本
docker-compose run --rm autosarm:20250120
```

## 镜像管理

### 查看所有镜像

```bash
# 查看 autosarm 的所有标签
docker images autosarm

# 格式化输出
docker images autosarm --format "{{.Repository}}:{{.Tag}} ({{.Size}}, {{.CreatedAt}})"
```

### 清理旧镜像

```bash
# 删除特定日期的镜像
docker rmi autosarm:20250120

# 只保留最近 7 天的镜像
docker images autosarm --format "{{.Tag}}" | \
  grep -E '^[0-9]{8}$' | \
  sort -r | \
  tail -n +8 | \
  xargs -I {} docker rmi autosarm:{}

# 清理未使用的镜像
docker image prune -f
```

### 导出和导入镜像

```bash
# 导出特定日期的镜像
docker save autosarm:20250124 | gzip > autosarm-20250124.tar.gz

# 导入镜像
gunzip -c autosarm-20250124.tar.gz | docker load
```

## 最佳实践

### 1. 生产环境使用

在生产环境中，建议使用日期标签而不是 `latest`：

```yaml
# docker-compose.prod.yml
services:
  autosarm:
    image: autosarm:20250124  # 使用固定日期标签
    # ... 其他配置
```

### 2. 开发环境使用

在开发环境中，可以使用 `latest` 标签：

```yaml
# docker-compose.yml
services:
  autosarm:
    image: autosarm:latest
    # ... 其他配置
```

### 3. 版本追踪

在项目文档中记录关键版本：

```markdown
## 版本历史

- 20250124: 添加新功能 X
- 20250120: 修复 bug Y
- 20250115: 初始版本
```

### 4. 自动化构建

在 CI/CD 流程中使用：

```bash
#!/bin/bash
DATE_TAG=$(date +%Y%m%d)
BUILD_NUMBER=$CI_BUILD_NUMBER

# 构建并推送
docker build -t myregistry/autosarm:latest \
             -t myregistry/autosarm:$DATE_TAG \
             -t myregistry/autosarm:build-$BUILD_NUMBER \
             .

docker push myregistry/autosarm:latest
docker push myregistry/autosarm:$DATE_TAG
```

## 常见问题

### Q: 如何知道 latest 对应哪个日期？

```bash
# 查看镜像 ID
docker images autosarm:latest --format "{{.ID}}"

# 查找相同 ID 的其他标签
docker images autosarm --format "{{.ID}} {{.Tag}}" | grep <IMAGE_ID>
```

### Q: 如何在两个版本之间对比？

```bash
# 对比镜像历史
docker history autosarm:20250124
docker history autosarm:20250120

# 运行并对比
docker run --rm autosarm:20250124 pip list > packages-20250124.txt
docker run --rm autosarm:20250120 pip list > packages-20250120.txt
diff packages-20250124.txt packages-20250120.txt
```

### Q: 标签太多了怎么办？

定期清理旧标签，只保留必要的版本：

```bash
# 创建清理脚本
cat > cleanup_old_images.sh << 'EOF'
#!/bin/bash
KEEP_DAYS=30
CUTOFF_DATE=$(date -d "$KEEP_DAYS days ago" +%Y%m%d)

docker images autosarm --format "{{.Tag}}" | \
  grep -E '^[0-9]{8}$' | \
  while read tag; do
    if [ "$tag" -lt "$CUTOFF_DATE" ]; then
      echo "删除旧镜像: autosarm:$tag"
      docker rmi autosarm:$tag
    fi
  done
EOF

chmod +x cleanup_old_images.sh
./cleanup_old_images.sh
```

## 示例工作流

### 每日构建工作流

```bash
# 1. 构建新镜像
./quick_build.sh

# 2. 测试新镜像
./build.sh test

# 3. 如果测试通过，使用新镜像
docker-compose down
docker-compose up -d

# 4. 如果测试失败，回退到昨天的版本
YESTERDAY=$(date -d "yesterday" +%Y%m%d)
docker tag autosarm:$YESTERDAY autosarm:latest
docker-compose up -d
```

### 发布工作流

```bash
# 1. 构建发布版本
DATE_TAG=$(date +%Y%m%d)
docker build -t autosarm:latest \
             -t autosarm:$DATE_TAG \
             -t autosarm:v1.0.0 \
             .

# 2. 测试
./build.sh test

# 3. 推送到仓库（如果有）
docker push autosarm:$DATE_TAG
docker push autosarm:v1.0.0
docker push autosarm:latest

# 4. 创建 Git 标签
git tag -a v1.0.0-$DATE_TAG -m "Release v1.0.0 built on $DATE_TAG"
git push origin v1.0.0-$DATE_TAG
```
