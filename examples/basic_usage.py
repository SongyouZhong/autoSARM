#!/usr/bin/env python
"""
Basic AutoSARM Usage Example

This script demonstrates the basic workflow for generating SAR matrices
from a compound library.
"""

import pandas as pd
from pathlib import Path

# Import from the autosarm package
from autosarm import fragmentize, create_sarm_matrix
from autosarm.utils import canonic_smiles, df_valid, float_row


def main():
    """Run basic SAR matrix generation example."""
    
    print("=" * 60)
    print("AutoSARM Basic Usage Example")
    print("=" * 60)
    
    # Example compound data
    # In real usage, you would load from a CSV file
    data = {
        'smiles': [
            'Cc1ccc(C2CCN(C(=O)c3ccccc3)CC2)cc1',
            'Cc1ccc(C2CCN(C(=O)c3cccc(F)c3)CC2)cc1',
            'Cc1ccc(C2CCN(C(=O)c3cccc(Cl)c3)CC2)cc1',
            'CCc1ccc(C2CCN(C(=O)c3ccccc3)CC2)cc1',
            'Fc1ccc(C2CCN(C(=O)c3ccccc3)CC2)cc1',
            'Cc1ccc(C2CCN(C(=O)c3ccc(F)cc3)CC2)cc1',
            'Cc1ccc(C2CCN(C(=O)c3ccc(Cl)cc3)CC2)cc1',
            'Cc1ccc(C2CCN(C(=O)c3ccncc3)CC2)cc1',
        ],
        'IC50_uM': [0.15, 0.08, 0.12, 0.45, 0.22, 0.19, 0.14, 0.89],
        'compound_id': [f'COMP{i:03d}' for i in range(1, 9)]
    }
    
    df = pd.DataFrame(data)
    
    print(f"\nLoaded {len(df)} compounds")
    print("\nSample data:")
    print(df.head())
    
    # Step 1: Validate SMILES
    print("\n" + "-" * 40)
    print("Step 1: Validating SMILES")
    df_clean = df_valid(df, row_smi='smiles')
    print(f"Valid compounds: {len(df_clean)}")
    
    # Step 2: Canonicalize SMILES
    print("\n" + "-" * 40)
    print("Step 2: Canonicalizing SMILES")
    df_clean['Cano_SMILES'] = df_clean['smiles'].apply(canonic_smiles)
    df_clean = df_clean.drop_duplicates(subset=['Cano_SMILES'])
    print(f"Unique compounds: {len(df_clean)}")
    
    # Step 3: Convert activity to float
    print("\n" + "-" * 40)
    print("Step 3: Processing activity data")
    df_clean = float_row(df_clean, cols=['IC50_uM'])
    
    # Step 4: Fragmentize molecules
    print("\n" + "-" * 40)
    print("Step 4: Fragmentizing molecules")
    round1, round2 = fragmentize(
        df_clean['Cano_SMILES'],
        n_jobs=1,  # Use 1 for this example
        drop_duplicate=False,
        pos_args={'RR': True, 'nRnR': True}
    )
    print(f"Round 1 fragments: {len(round1)}")
    print(f"Round 2 fragments: {len(round2)}")
    
    # Step 5: Create SAR matrices
    print("\n" + "-" * 40)
    print("Step 5: Creating SAR matrices")
    
    output_folder = Path("./SAR_Example_Results")
    output_folder.mkdir(exist_ok=True)
    
    left_info, right_info, combined_info = create_sarm_matrix(
        round1, round2, df_clean,
        save_folder=str(output_folder),
        smi_col='Cano_SMILES',
        value_col=['IC50_uM'],
        minimum_count_site1=2,
        minimum_count_site2=2,
        csv2excel=False,
        n_jobs=1
    )
    
    # Results summary
    print("\n" + "=" * 60)
    print("Results Summary")
    print("=" * 60)
    print(f"Left tables created: {len(left_info)}")
    print(f"Right tables created: {len(right_info)}")
    print(f"Combined tables created: {len(combined_info)}")
    print(f"\nResults saved to: {output_folder.absolute()}")
    
    # List output files
    print("\nOutput files:")
    for f in sorted(output_folder.rglob("*.csv")):
        print(f"  - {f.relative_to(output_folder)}")
    
    print("\n" + "=" * 60)
    print("Example complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
