# AutoSARM

**Automatic Structure-Activity Relationship Matrix Generator**

A comprehensive toolkit for drug discovery research that enables automated SAR analysis and molecular fragmentation.

[![Python](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/)
[![RDKit](https://img.shields.io/badge/RDKit-Required-green.svg)](https://www.rdkit.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [API Reference](#api-reference)
- [Contributing](#contributing)
- [License](#license)

---

## 🎯 Overview

AutoSARM is a professional computational chemistry tool designed to help medicinal chemists and researchers rapidly analyze and visualize structure-activity relationships in compound series. Through automated molecular fragmentation and statistical analysis, this tool enables:

- 🔬 Automatic identification of molecular core scaffolds and substituents
- 📊 Generation of 1D and 2D SAR tables
- 🌳 Construction of molecular derivatization trees
- 📈 Visualization of SAR patterns
- 💊 Support for lead compound optimization

---

## ✨ Features

### Core Capabilities

1. **Intelligent Molecular Fragmentation**
   - Ring-adjacent single bond cleavage
   - Ring-ring bond cleavage (optional)
   - Two-round iterative fragmentation

2. **SAR Matrix Generation**
   - Left/Right attachment site analysis
   - Combined multi-site analysis
   - Single-cut fragment analysis

3. **SAR Tree Visualization**
   - Hierarchical structure display
   - Activity-based highlighting
   - Graphviz-powered visualization

4. **Activity Statistics**
   - Mean, median, std calculations
   - Min/max range analysis
   - Substructure-based matching

---

## 📦 Installation

### From Source (Recommended)

```bash
# Clone the repository
git clone https://github.com/yourorg/autosarm.git
cd autosarm

# Create conda environment
conda env create -f environment.yml
conda activate autosarm

# Install in development mode
pip install -e .
```

### Using pip

```bash
pip install autosarm
```

### Dependencies

Core dependencies:
- Python >= 3.9
- RDKit >= 2023.03.1
- pandas >= 1.5.0
- numpy >= 1.21.0
- pandarallel >= 1.6.0
- graphviz >= 0.20.0

---

## 🚀 Quick Start

### Command Line Interface

```bash
# Create SAR matrix
autosarm sarm --csvFile compounds.csv --column IC50_uM --save_folder SAR_Results

# Create SAR tree
autosarm tree --fragment_core "c1ccc(*)*cc1" --workFolder SAR_Results --rootTitle MyTree

# Get help
autosarm --help
```

### Python API

```python
from autosarm import fragmentize, create_sarm_matrix
from autosarm.utils import canonic_smiles
import pandas as pd

# Load your data
df = pd.read_csv("compounds.csv")
df['Cano_SMILES'] = df['smiles'].apply(canonic_smiles)

# Fragmentize molecules
round1, round2 = fragmentize(df['Cano_SMILES'], n_jobs=8)

# Create SAR matrices
left, right, combined = create_sarm_matrix(
    round1, round2, df,
    save_folder='SAR_Results',
    smi_col='Cano_SMILES',
    value_col=['IC50_uM']
)
```

---

## 📖 Usage

### Input Data Format

Your input CSV file should contain:

| Column | Required | Description |
|--------|----------|-------------|
| smiles | Yes | Valid SMILES strings |
| activity_column | Yes | Numeric activity values |
| compound_id | No | Compound identifiers |

Example:
```csv
smiles,IC50_uM,compound_id
Cc1ccc(C2CCN(C(=O)c3ccccc3)CC2)cc1,0.15,COMP001
Cc1ccc(C2CCN(C(=O)c3cccc(F)c3)CC2)cc1,0.08,COMP002
```

### Command Line Options

#### `autosarm sarm`

| Option | Default | Description |
|--------|---------|-------------|
| `--csvFile` | Required | Input CSV file path |
| `--column` | Required | Activity column name(s) |
| `--type` | smiles | Analysis type (smiles/scaffold) |
| `--minimumSite1` | 3 | Min fragment count for site 1 |
| `--minimumSite2` | 3 | Min fragment count for site 2 |
| `--n_jobs` | 8 | Number of parallel jobs |
| `--save_folder` | SAR_Results | Output directory |
| `--csv2excel` | 0 | Generate Excel files (0/1) |

#### `autosarm tree`

| Option | Default | Description |
|--------|---------|-------------|
| `--fragment_core` | Required | Core fragment SMILES |
| `--rootTitle` | Required | Tree title |
| `--workFolder` | Required | Working directory |
| `--maxLevel` | 5 | Maximum tree depth |
| `--treeContent` | ['double-cut'] | Tree content types |

---

## 📁 Project Structure

```
autosarm/
├── pyproject.toml          # Package configuration
├── environment.yml         # Conda environment
├── README.md              # This file
├── src/
│   └── autosarm/          # Main package
│       ├── __init__.py
│       ├── core/          # Core algorithms
│       │   ├── sarm.py    # SAR matrix generation
│       │   └── tree.py    # SAR tree generation
│       ├── cli/           # Command line interface
│       │   ├── main.py
│       │   ├── create_sarm.py
│       │   └── create_tree.py
│       └── utils/         # Utility functions
│           ├── mol_utils.py
│           └── data_utils.py
├── tests/                 # Test suite
│   ├── test_mol_utils.py
│   └── test_sarm.py
├── examples/              # Example scripts
│   ├── basic_usage.py
│   └── example_compounds.csv
├── docs/                  # Documentation
└── scripts/               # Utility scripts
```

---

## 📚 API Reference

### Core Functions

#### `fragmentize(smiles_list, n_jobs=20, pos_args=None)`
Perform two-round molecular fragmentation.

#### `create_sarm_matrix(round1, round2, df_active, save_folder, ...)`
Generate SAR matrices from fragmentation data.

#### `create_sar_tree(fragment_core, root_title, work_folder, ...)`
Create SAR tree visualization.

### Utility Functions

#### `get_mol(smiles_or_mol)`
Convert SMILES to RDKit Mol object.

#### `canonic_smiles(smiles_or_mol)`
Canonicalize SMILES string.

#### `compute_fingerprint(mol, radius=2, n_bits=1024)`
Compute Morgan fingerprint.

#### `compute_similarity(smi1, smi2, mode='smi-smi')`
Calculate Tanimoto similarity.

---

## 🧪 Running Tests

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest tests/ --cov=autosarm --cov-report=html

# Run specific test file
pytest tests/test_mol_utils.py -v
```

---

## 🐳 Docker Usage

```bash
# Build image
docker build -t autosarm .

# Run container
docker run -v $(pwd)/data:/app/data autosarm \
    autosarm sarm --csvFile /app/data/compounds.csv --column IC50_uM
```

See [DOCKER_USAGE.md](DOCKER_USAGE.md) for detailed Docker instructions.

---

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## 📄 License

This project is licensed under the MIT License (declared in `pyproject.toml`). A standalone `LICENSE` file has not yet been added to the repo.

---

## 📧 Contact

For questions or support, please open an issue on GitHub.

---

## 🙏 Acknowledgments

- RDKit development team
- All contributors to this project
