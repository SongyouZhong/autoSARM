# AutoSARM Examples

This directory contains example scripts and notebooks demonstrating how to use AutoSARM.

## Contents

- `basic_usage.py` - Basic example of SAR matrix generation
- `example_compounds.csv` - Sample compound dataset
- `advanced_usage.ipynb` - Jupyter notebook with advanced features (optional)

## Quick Start

1. Install the package:
```bash
pip install -e ..
```

2. Run the basic example:
```bash
python basic_usage.py
```

## Example Workflow

### 1. Prepare Your Data

Your input CSV should contain at minimum:
- `smiles` column with valid SMILES strings
- One or more activity columns (e.g., `IC50_uM`)

### 2. Generate SAR Matrix

```python
from autosarm import create_sarm_matrix, fragmentize
from autosarm.utils import canonic_smiles, df_valid
import pandas as pd

# Load data
df = pd.read_csv("your_compounds.csv")

# Validate and prepare
df_valid = df_valid(df, row_smi='smiles')
df['Cano_SMILES'] = df['smiles'].apply(canonic_smiles)

# Fragmentize
round1, round2 = fragmentize(df['Cano_SMILES'], n_jobs=4)

# Create matrices
left, right, combined = create_sarm_matrix(
    round1, round2, df,
    save_folder='output',
    smi_col='Cano_SMILES',
    value_col=['IC50_uM']
)
```

### 3. Generate SAR Tree

```python
from autosarm import create_sar_tree

create_sar_tree(
    fragment_core="c1ccc(*)*cc1",  # Your core structure
    root_title="Example_Tree",
    work_folder="output",
    max_level=5
)
```

## Tips

- Use `n_jobs` parameter to control parallelization
- Set `minimum_count_site1/2` to filter rare fragments
- Enable `csv2excel` for visual Excel reports
