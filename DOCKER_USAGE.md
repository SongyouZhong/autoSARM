# AutoSARM Docker Usage Guide

This guide explains how to use AutoSARM with Docker.

## Prerequisites

- Docker >= 20.10
- Docker Compose >= 2.0

## Quick Start

### 1. Build the Image

```bash
# Build production image
./docker-manage.sh build

# Or build development image
./docker-manage.sh build development
```

### 2. Run SAR Analysis

```bash
# Place your data in the data/ directory
cp your_compounds.csv data/

# Run SAR matrix generation
./docker-manage.sh run sarm --csvFile /app/data/your_compounds.csv --column IC50_uM --save_folder /app/output
```

### 3. View Results

Results will be saved to the `output/` directory.

## Commands

### Using docker-manage.sh

```bash
# Build image
./docker-manage.sh build [production|development]

# Run SAR generation
./docker-manage.sh run sarm --csvFile /app/data/compounds.csv --column IC50_uM

# Start development shell
./docker-manage.sh dev

# Run tests
./docker-manage.sh test

# Start Jupyter notebook
./docker-manage.sh jupyter

# Show status
./docker-manage.sh status

# Clean up
./docker-manage.sh clean
```

### Using Docker Compose Directly

```bash
# Build
docker compose build

# Run SAR generation
docker compose run --rm autosarm sarm \
    --csvFile /app/data/compounds.csv \
    --column IC50_uM \
    --save_folder /app/output

# Start development container
docker compose run --rm autosarm-dev

# Run tests
docker compose run --rm autosarm-test

# Start Jupyter
docker compose up autosarm-jupyter
```

## Directory Structure

```
autoSARM/
├── data/           # Input data (mounted to /app/data)
├── output/         # Output results (mounted to /app/output)
├── notebooks/      # Jupyter notebooks (mounted to /app/notebooks)
└── ...
```

## Examples

### Example 1: Basic SAR Matrix Generation

```bash
# Prepare data
mkdir -p data output
cp examples/example_compounds.csv data/

# Run analysis
./docker-manage.sh run sarm \
    --csvFile /app/data/example_compounds.csv \
    --column IC50_uM \
    --save_folder /app/output/results \
    --minimumSite1 2 \
    --minimumSite2 2

# Check results
ls -la output/results/
```

### Example 2: Create SAR Tree

```bash
# First run SAR matrix generation, then create tree
./docker-manage.sh run tree \
    --fragment_core "c1ccc(*)cc1*" \
    --workFolder /app/output/results \
    --rootTitle MyTree \
    --maxLevel 5
```

### Example 3: Interactive Development

```bash
# Start development shell
./docker-manage.sh dev

# Inside container:
python -c "from autosarm import fragmentize; print('OK')"
pytest tests/ -v
python examples/basic_usage.py
```

### Example 4: Jupyter Notebook

```bash
# Start Jupyter server
./docker-manage.sh jupyter

# Open browser at http://localhost:8888
# Create notebooks in the notebooks/ directory
```

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGE_TAG` | `latest` | Docker image tag |
| `JUPYTER_PORT` | `8888` | Jupyter notebook port |

Example:
```bash
IMAGE_TAG=v1.0 ./docker-manage.sh build
JUPYTER_PORT=9999 ./docker-manage.sh jupyter
```

## Troubleshooting

### Permission Issues

If you encounter permission issues with output files:

```bash
# Fix ownership
sudo chown -R $USER:$USER output/
```

### Build Failures

```bash
# Clean and rebuild
./docker-manage.sh clean
./docker-manage.sh build
```

### Memory Issues

For large datasets, increase Docker memory limit:

```bash
# In Docker Desktop settings, increase memory allocation
# Or use command line:
docker compose run --rm -m 8g autosarm sarm --csvFile ...
```

## Multi-Platform Build

To build for multiple platforms:

```bash
docker buildx build --platform linux/amd64,linux/arm64 -t autosarm:latest .
```

## CI/CD Integration

Example GitHub Actions workflow:

```yaml
name: Build and Test

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      
      - name: Build Docker image
        run: docker compose build autosarm-test
      
      - name: Run tests
        run: docker compose run --rm autosarm-test
```
