#!/bin/bash
# =============================================================================
# AutoSARM Docker Management Script
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
    echo "  build [target]    Build Docker image (targets: production, development)"
    echo "  run [args]        Run SAR matrix generation with arguments"
    echo "  dev               Start development container with shell"
    echo "  test              Run test suite"
    echo "  jupyter           Start Jupyter notebook server"
    echo "  shell             Start interactive shell in production container"
    echo "  clean             Remove containers and images"
    echo "  logs [service]    Show logs for a service"
    echo "  status            Show status of containers"
    echo ""
    echo "Examples:"
    echo "  $0 build                                    # Build production image"
    echo "  $0 build development                        # Build development image"
    echo "  $0 run sarm --csvFile data/compounds.csv --column IC50_uM"
    echo "  $0 dev                                      # Start dev shell"
    echo "  $0 test                                     # Run tests"
    echo "  $0 jupyter                                  # Start Jupyter"
    echo ""
}

# Build Docker image
cmd_build() {
    local target="${1:-production}"
    print_msg $BLUE "Building AutoSARM Docker image (target: $target)..."
    
    docker compose build --build-arg BUILDKIT_INLINE_CACHE=1 \
        $([ "$target" = "development" ] && echo "autosarm-dev" || echo "autosarm")
    
    print_msg $GREEN "Build complete!"
}

# Run SAR generation
cmd_run() {
    print_msg $BLUE "Running AutoSARM..."
    
    # Create output directory if not exists
    mkdir -p output data
    
    docker compose run --rm autosarm "$@"
}

# Start development container
cmd_dev() {
    print_msg $BLUE "Starting development container..."
    
    mkdir -p output data notebooks
    
    docker compose run --rm autosarm-dev
}

# Run tests
cmd_test() {
    print_msg $BLUE "Running tests..."
    
    docker compose run --rm autosarm-test
    
    print_msg $GREEN "Tests complete!"
}

# Start Jupyter
cmd_jupyter() {
    print_msg $BLUE "Starting Jupyter Notebook..."
    print_msg $YELLOW "Access at: http://localhost:${JUPYTER_PORT:-8888}"
    
    mkdir -p output data notebooks
    
    docker compose up autosarm-jupyter
}

# Start shell
cmd_shell() {
    print_msg $BLUE "Starting interactive shell..."
    
    docker compose run --rm --entrypoint /bin/bash autosarm
}

# Clean up
cmd_clean() {
    print_msg $YELLOW "Cleaning up Docker resources..."
    
    docker compose down --rmi local --volumes --remove-orphans 2>/dev/null || true
    docker image prune -f
    
    print_msg $GREEN "Cleanup complete!"
}

# Show logs
cmd_logs() {
    local service="${1:-autosarm}"
    docker compose logs -f "$service"
}

# Show status
cmd_status() {
    print_msg $BLUE "Container Status:"
    docker compose ps
    
    echo ""
    print_msg $BLUE "Images:"
    docker images | grep autosarm || echo "No autosarm images found"
}

# Main
main() {
    local cmd="${1:-}"
    shift || true
    
    case "$cmd" in
        build)
            cmd_build "$@"
            ;;
        run)
            cmd_run "$@"
            ;;
        dev)
            cmd_dev
            ;;
        test)
            cmd_test
            ;;
        jupyter)
            cmd_jupyter
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
        status)
            cmd_status
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
