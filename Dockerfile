# 使用 micromamba 基础镜像
FROM mambaorg/micromamba:1.5.8

# 切换到 root 用户以安装系统包
USER root

# 修复 /tmp 目录权限
RUN chmod 1777 /tmp

# 安装系统依赖（graphviz 需要系统级安装）
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    graphviz \
    libgraphviz-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# 设置工作目录
WORKDIR /app

# 复制环境配置文件
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml

# 使用 micromamba 创建环境
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# 切换到非 root 用户
USER $MAMBA_USER

# 复制项目文件
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

# 激活环境并设置环境变量
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# 设置 Python 路径
ENV PYTHONPATH=/app

# 暴露端口（如果需要 Jupyter）
EXPOSE 8888

# 默认命令（可以根据需要修改）
CMD ["/bin/bash"]
