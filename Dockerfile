# =============================================================================
# AutoSARM Docker Image
# Multi-stage build for efficient production deployment
# =============================================================================

# Stage 1: Build stage
FROM mambaorg/micromamba:1.5.8 AS builder

USER root

# Fix /tmp permissions
RUN chmod 1777 /tmp

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    graphviz \
    libgraphviz-dev \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy environment file first (for better caching)
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Create conda environment
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Copy project files
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml /app/
COPY --chown=$MAMBA_USER:$MAMBA_USER src/ /app/src/

# Install the package in development mode
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN micromamba run -n base pip install -e . --no-deps

# =============================================================================
# Stage 2: Production stage
# =============================================================================
FROM mambaorg/micromamba:1.5.8 AS production

USER root

# Fix /tmp permissions
RUN chmod 1777 /tmp

# Install runtime dependencies (graphviz for rendering, curl for healthcheck)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    graphviz \
    curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy environment from builder
COPY --from=builder /opt/conda /opt/conda

WORKDIR /app

# Copy installed package
COPY --from=builder --chown=$MAMBA_USER:$MAMBA_USER /app /app

# Copy additional files
COPY --chown=$MAMBA_USER:$MAMBA_USER README.md /app/
COPY --chown=$MAMBA_USER:$MAMBA_USER examples/ /app/examples/

# Copy docker entrypoint
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/docker-entrypoint.sh /app/docker-entrypoint.sh
RUN chmod +x /app/docker-entrypoint.sh

# Create data directories (before switching to non-root user)
RUN mkdir -p /app/data /app/output /app/logs /tmp/autosarm && \
    chown -R $MAMBA_USER:$MAMBA_USER /app/data /app/output /app/logs /tmp/autosarm

# Switch to non-root user
USER $MAMBA_USER

# Activate environment
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Set Python path
ENV PYTHONPATH=/app/src

# Expose API port
EXPOSE 8030

# Health check — hit the /health endpoint
HEALTHCHECK --interval=30s --timeout=10s --start-period=15s --retries=3 \
    CMD curl -f http://localhost:${PORT:-8030}/health || exit 1

# Default: start worker (API + DB polling)
ENTRYPOINT ["/app/docker-entrypoint.sh"]
CMD ["serve"]
