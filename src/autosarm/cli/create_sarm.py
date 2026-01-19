"""
CLI module for creating SAR matrices.
"""

from __future__ import annotations

import argparse
import logging
import math
import os
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

from autosarm.core.sarm import create_sarm_matrix, fragmentize
from autosarm.utils.mol_utils import canonic_smiles
from autosarm.utils.data_utils import float_row, df_valid


def setup_logging(save_folder: str) -> str:
    """
    Set up logging configuration.
    
    Args:
        save_folder: Directory for log files
        
    Returns:
        Path to log file
    """
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(save_folder, f'create_sarm_{timestamp}.log')
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, encoding='utf-8'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    logging.info(f"Log file created: {log_file}")
    return log_file


def run_sarm(args) -> None:
    """
    Run SAR matrix creation.
    
    Args:
        args: Parsed command line arguments
    """
    try:
        # Setup logging
        log_file = setup_logging(args.save_folder)
        logging.info("=" * 80)
        logging.info("Starting SAR Matrix Generation")
        logging.info("=" * 80)
        
        act_csv_file = args.csvFile
        value_cols = args.column
        log_transform = args.log
        
        logging.info(f"Input CSV file: {act_csv_file}")
        logging.info(f"Value columns: {value_cols}")
        logging.info(f"Log transform: {log_transform}")
        logging.info(f"Minimum count site1: {args.minimumSite1}")
        logging.info(f"Minimum count site2: {args.minimumSite2}")
        logging.info(f"Number of jobs: {args.n_jobs}")
        logging.info(f"Save folder: {args.save_folder}")
        
        # Load and process data
        df_act = pd.read_csv(act_csv_file)
        logging.info(f"Loaded {len(df_act)} rows from CSV file")
        
        df_act = float_row(df_act, cols=value_cols)
        logging.info("Converted value columns to float")
        
        if log_transform:
            for col in value_cols:
                df_act[col] = df_act[col].apply(lambda x: math.log(x) if x > 0 else x)
            logging.info("Applied log transform")
        
        df_sele = df_valid(df_act, row_smi='smiles')
        logging.info(f"Validated SMILES, {len(df_sele)} valid molecules")
        
        if 'Cano_SMILES' not in df_sele.columns:
            logging.info("Canonicalizing SMILES...")
            df_sele['Cano_SMILES'] = df_sele['smiles'].apply(canonic_smiles)
            logging.info("SMILES canonicalization complete")
        
        df_sele = df_sele.drop_duplicates(subset=['Cano_SMILES'])
        logging.info(f"After removing duplicates: {len(df_sele)} unique molecules")
        
        if args.type == 'scaffold':
            from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
            logging.info("Extracting Murcko scaffolds...")
            df_sele['Scaffold'] = [MurckoScaffoldSmiles(smi) for smi in df_sele["Cano_SMILES"]]
            act_cpds = df_sele['Scaffold']
            smi_col = 'Scaffold'
            logging.info(f"Extracted {len(df_sele['Scaffold'].unique())} unique scaffolds")
        else:
            act_cpds = df_sele['Cano_SMILES']
            smi_col = 'Cano_SMILES'
            logging.info("Using full SMILES for analysis")
        
        logging.info("Starting fragmentization...")
        df_round1, df_round2 = fragmentize(
            act_cpds, 
            n_jobs=args.n_jobs, 
            drop_duplicate=False, 
            pos_args={'RR': True, 'nRnR': True}
        )
        logging.info(f"Round 1 fragments: {len(df_round1)}")
        logging.info(f"Round 2 fragments: {len(df_round2)}")
        
        logging.info("Creating SARM tables...")
        df_left, df_right, df_combine = create_sarm_matrix(
            df_round1, df_round2, df_sele,
            save_folder=args.save_folder,
            smi_col=smi_col,
            value_col=value_cols,
            minimum_count_site1=int(args.minimumSite1),
            minimum_count_site2=int(args.minimumSite2),
            csv2excel=bool(args.csv2excel),
            cal_table_stats=True,
            n_jobs=args.n_jobs
        )
        
        logging.info("=" * 80)
        logging.info("SAR Matrix Generation Completed Successfully!")
        logging.info(f"Left tables: {len(df_left)}")
        logging.info(f"Right tables: {len(df_right)}")
        logging.info(f"Combined tables: {len(df_combine)}")
        logging.info(f"Results saved to: {args.save_folder}")
        logging.info(f"Log file: {log_file}")
        logging.info("=" * 80)
        
    except Exception as e:
        logging.error("=" * 80)
        logging.error("ERROR OCCURRED!")
        logging.error("=" * 80)
        logging.error(f"Error type: {type(e).__name__}")
        logging.error(f"Error message: {str(e)}")
        logging.error("Full traceback:", exc_info=True)
        logging.error("=" * 80)
        raise


def get_parser() -> argparse.ArgumentParser:
    """Create argument parser for standalone usage."""
    parser = argparse.ArgumentParser(
        description="Create SAR matrix from compound library"
    )
    parser.add_argument("--csvFile", required=True, help="Input CSV file with SMILES")
    parser.add_argument("--type", default="smiles", choices=["smiles", "scaffold"])
    parser.add_argument("--column", nargs="+", required=True, help="Activity column(s)")
    parser.add_argument("--log", type=int, default=0)
    parser.add_argument("--minimumSite1", type=float, default=3)
    parser.add_argument("--minimumSite2", type=float, default=3)
    parser.add_argument("--n_jobs", type=int, default=8)
    parser.add_argument("--save_folder", default="SAR_Results")
    parser.add_argument("--csv2excel", type=int, default=0)
    return parser


def main():
    """Main entry point for standalone script."""
    parser = get_parser()
    args = parser.parse_args()
    run_sarm(args)


if __name__ == "__main__":
    main()
