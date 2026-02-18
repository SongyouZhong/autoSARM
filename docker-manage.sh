#!/bin/bash
# =============================================================================
# AutoSARM Docker Management Script
# 支持多实例 Worker 模式
# =============================================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values
IMAGE_TAG="${IMAGE_TAG:-latest}"
COMPOSE_FILE="docker-compose.yml"

# Print colored message
print_msg() {
    local color=$1
    local msg=$2
    echo -e "${color}${msg}${NC}"
}

# Print usage
usage() {
    echo "AutoSARM Docker Management Script"
    echo ""
    echo "Usage: $0 <command> [options]"
    echo ""
    echo "Commands:"
    echo "  build             Build Docker image"
    echo "  start [n]         Start n worker instances (default: all 3)"
    echo "  stop              Stop all workers"
    echo "  restart           Restart all workers"
    echo "  run [args]        Run SAR matrix generation locally (CLI mode)"
    echo "  shell             Start interactive shell in a worker container"
    echo "  clean             Remove containers and images"
    echo "  logs [service]    Show logs (default: all workers)"
    echo "  status            Show status of containers"
    echo "  health            Check health of all workers"
    echo ""
    echo "Examples:"
    echo "  $0 build                                    # Build image"
    echo "  $0 start                                    # Start all 3 workers"
    echo "  $0 start 1                                  # Start only worker 1"
    echo "  $0 start 2                                  # Start workers 1 and 2"
    echo "  $0 stop                                     # Stop all workers"
    echo "  $0 logs autosarm-1                          # Logs for worker 1"
    echo "  $0 health                                   # Check health endpoints"
    echo "  $0 run sarm --csvFile data/compounds.csv --column IC50_uM"
    echo ""
}

# Build Docker image
cmd_build() {
    print_msg $BLUE "Building AutoSARM Docker image..."

    docker compose -f "$COMPOSE_FILE" build

    print_msg $GREEN "Build complete!"
}

# Start workers
cmd_start() {
    local count="${1:-3}"
    print_msg $BLUE "Starting $count AutoSARM worker(s)..."

    mkdir -p output data

    case "$count" in
        1)
            docker compose -f "$COMPOSE_FILE" up -d autosarm-1
            ;;
        2)
            docker compose -f "$COMPOSE_FILE" up -d autosarm-1 autosarm-2
            ;;
        *)
            docker compose -f "$COMPOSE_FILE" up -d
            ;;
    esac

    print_msg $GREEN "Workers started!"
    echo ""
    cmd_status
}

# Stop workers
cmd_stop() {
    print_msg $YELLOW "Stopping all workers..."
    docker compose -f "$COMPOSE_FILE" down
    print_msg $GREEN "All workers stopped"
}

# Restart workers
cmd_restart() {
    cmd_stop
    sleep 2
    cmd_start "${1:-3}"
}

# Run SAR generation (CLI mode, one-off)
cmd_run() {
    print_msg $BLUE "Running AutoSARM (CLI mode)..."

    mkdir -p output data

    docker compose -f "$COMPOSE_FILE" run --rm autosarm-1 "$@"
}

# Start shell
cmd_shell() {
    print_msg $BLUE "Starting interactive shell..."

    docker compose -f "$COMPOSE_FILE" run --rm autosarm-1 shell
}

# Clean up
cmd_clean() {
    print_msg $YELLOW "Cleaning up Docker resources..."

    docker compose -f "$COMPOSE_FILE" down --rmi local --volumes --remove-orphans 2>/dev/null || true
    docker image prune -f

    print_msg $GREEN "Cleanup complete!"
}

# Show logs
cmd_logs() {
    local service="${1:-}"
    if [ -n "$service" ]; then
        docker compose -f "$COMPOSE_FILE" logs -f "$service"
    else
        docker compose -f "$COMPOSE_FILE" logs -f
    fi
}

# Show status
cmd_status() {
    print_msg $BLUE "Container Status:"
    docker compose -f "$COMPOSE_FILE" ps

    echo ""
    print_msg $BLUE "Images:"
    docker images | grep autosarm || echo "No autosarm images found"
}

# Health check
cmd_health() {
    print_msg $BLUE "Checking worker health..."

    local ports=("${SARM_PORT_1:-8030}" "${SARM_PORT_2:-8031}" "${SARM_PORT_3:-8032}")
    local names=("autosarm-1" "autosarm-2" "autosarm-3")

    for i in 0 1 2; do
        local port="${ports[$i]}"
        local name="${names[$i]}"
        local response

        response=$(curl -s -o /dev/null -w "%{http_code}" "http://localhost:${port}/health" 2>/dev/null) || response="000"

        if [ "$response" = "200" ]; then
            print_msg $GREEN "  $name (port $port): healthy"
        else
            print_msg $RED "  $name (port $port): unreachable (HTTP $response)"
        fi
    done
}

# Main
main() {
    local cmd="${1:-}"
    shift || true

    case "$cmd" in
        build)
            cmd_build
            ;;
        start|up)
            cmd_start "$@"
            ;;
        stop|down)
            cmd_stop
            ;;
        restart)
            cmd_restart "$@"
            ;;
        run)
            cmd_run "$@"
            ;;
        shell)
            cmd_shell
            ;;
        clean)
            cmd_clean
            ;;
        logs)
            cmd_logs "$@"
            ;;
        status|ps)
            cmd_status
            ;;
        health)
            cmd_health
            ;;
        help|--help|-h)
            usage
            ;;
        *)
            print_msg $RED "Error: Unknown command '$cmd'"
            echo ""
            usage
            exit 1
            ;;
    esac
}

main "$@"
