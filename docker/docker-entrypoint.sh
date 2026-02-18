#!/bin/bash
set -e

# Docker entrypoint script for AutoSARM Service

echo "=========================================="
echo "AutoSARM Service"
echo "=========================================="
echo "Starting at: $(date)"
echo ""

# Initialize conda
source /opt/conda/etc/profile.d/conda.sh 2>/dev/null || \
eval "$(micromamba shell hook --shell bash)"
micromamba activate base 2>/dev/null || conda activate base 2>/dev/null || true

echo "Using Python: $(which python)"
echo "Python version: $(python --version)"
echo ""

# Default command
CMD=${1:-serve}

case "$CMD" in
    serve)
        echo "Starting AutoSARM Worker (API + DB polling)..."
        exec python -m autosarm serve \
            --host ${HOST:-0.0.0.0} \
            --port ${PORT:-8030}
        ;;

    sarm)
        echo "Running SARM matrix generation..."
        shift
        exec autosarm sarm "$@"
        ;;

    tree)
        echo "Running SAR tree generation..."
        shift
        exec autosarm tree "$@"
        ;;

    shell)
        echo "Starting shell..."
        exec /bin/bash
        ;;

    *)
        echo "Unknown command: $CMD"
        echo "Available commands: serve, sarm, tree, shell"
        exit 1
        ;;
esac
