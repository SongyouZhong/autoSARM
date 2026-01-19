# AutoSARM

This project has been refactored to follow Python packaging best practices.

## New Structure

The code has been reorganized into a proper Python package with the following structure:

```
autosarm/
├── pyproject.toml          # Modern Python package configuration
├── environment.yml         # Conda environment
├── README.md              # Documentation
├── src/
│   └── autosarm/          # Main package (pip installable)
│       ├── __init__.py    # Package exports
│       ├── core/          # Core algorithms
│       │   ├── sarm.py    # SAR matrix generation
│       │   └── tree.py    # SAR tree generation
│       ├── cli/           # Command line interface
│       │   ├── main.py    # Main CLI entry point
│       │   ├── create_sarm.py
│       │   ├── create_tree.py
│       │   └── print_sar.py
│       └── utils/         # Utility functions
│           ├── mol_utils.py   # Molecular operations
│           └── data_utils.py  # Data processing
├── tests/                 # pytest test suite
├── examples/              # Usage examples
├── docs/                  # Documentation (preserved)
└── scripts/               # Utility scripts (preserved)
```

## Key Changes

1. **Package Structure**: Code moved to `src/autosarm/` layout (PEP 517/518 compliant)

2. **Modern Configuration**: Uses `pyproject.toml` instead of `setup.py`

3. **CLI Entry Points**: Install provides `autosarm` command with subcommands:
   - `autosarm sarm` - Create SAR matrices
   - `autosarm tree` - Create SAR trees
   - `autosarm print` - Print SAR tables

4. **Clean Imports**: 
   ```python
   from autosarm import fragmentize, create_sarm_matrix
   from autosarm.utils import canonic_smiles, get_mol
   ```

5. **Type Hints**: Functions have proper type annotations

6. **Testing**: pytest-based test suite in `tests/`

7. **Examples**: Working examples in `examples/`

## Installation

```bash
# Development install
pip install -e .

# Or with conda
conda env create -f environment.yml
conda activate autosarm
pip install -e .
```

## Migration from Old Code

Old imports:
```python
from utils.sarm_utils import create_SARM, fragmentize
from utils.common_utils import get_mol, canonic_smiles
```

New imports:
```python
from autosarm.core import create_sarm_matrix, fragmentize
from autosarm.utils import get_mol, canonic_smiles
```

## Legacy Files

The original files (`create_sarm.py`, `create_tree.py`, `print_sar_single_core.py`, `utils/`) are preserved for reference but are no longer the primary code.
