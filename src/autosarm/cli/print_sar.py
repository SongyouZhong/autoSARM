"""
CLI module for printing SAR tables.
"""

from __future__ import annotations

import argparse
import logging
import sys

import pandas as pd

from autosarm.utils.mol_utils import get_mol, safe_mol_from_smiles
from autosarm.utils.data_utils import float_row


def run_print(args) -> None:
    """
    Run SAR table printing.
    
    Args:
        args: Parsed command line arguments
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    
    logging.info("Generating SAR Table")
    logging.info(f"Core structure: {args.core}")
    logging.info(f"Activity file: {args.actFile}")
    
    # Load activity data
    df_act = pd.read_csv(args.actFile)
    df_act = float_row(df_act, cols=args.actCol, dropna=False)
    
    # TODO: Implement full SAR table printing logic
    # This is a placeholder for the print_sar_single_core functionality
    
    logging.info(f"Loaded {len(df_act)} compounds")
    logging.info(f"Output will be saved to: {args.output}")
    
    # For now, just filter compounds matching the core
    from rdkit import Chem
    
    core_mol = Chem.MolFromSmarts(args.core)
    if core_mol is None:
        logging.error(f"Invalid core SMARTS: {args.core}")
        return
    
    matching_compounds = []
    for idx, row in df_act.iterrows():
        mol = safe_mol_from_smiles(row[args.smiCol])
        if mol is not None and mol.HasSubstructMatch(core_mol):
            matching_compounds.append(idx)
    
    df_match = df_act.loc[matching_compounds]
    logging.info(f"Found {len(df_match)} compounds matching core")
    
    df_match.to_csv(args.output, index=False)
    logging.info(f"Results saved to: {args.output}")


def get_parser() -> argparse.ArgumentParser:
    """Create argument parser for standalone usage."""
    parser = argparse.ArgumentParser(
        description="Print SAR table for a specific core"
    )
    parser.add_argument("--core", required=True)
    parser.add_argument("--actFile", required=True)
    parser.add_argument("--smiCol", default="smiles")
    parser.add_argument("--actCol", nargs="+", default=["IC50"])
    parser.add_argument("--output", default="sar_table.csv")
    return parser


def main():
    """Main entry point for standalone script."""
    parser = get_parser()
    args = parser.parse_args()
    run_print(args)


if __name__ == "__main__":
    main()
