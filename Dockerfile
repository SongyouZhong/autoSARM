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

# Install runtime dependencies only
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    graphviz \
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

# Create data directories (before switching to non-root user)
RUN mkdir -p /app/data /app/output && \
    chown -R $MAMBA_USER:$MAMBA_USER /app/data /app/output

# Switch to non-root user
USER $MAMBA_USER

# Activate environment
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Set Python path
ENV PYTHONPATH=/app/src

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD micromamba run -n base python -c "import autosarm; print('OK')" || exit 1

# Default command - show help
ENTRYPOINT ["micromamba", "run", "-n", "base", "autosarm"]
CMD ["--help"]
