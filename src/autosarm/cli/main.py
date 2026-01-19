"""
Command Line Interface for AutoSARM.

This module provides the main CLI entry point for the autosarm package.
"""

from __future__ import annotations

import argparse
import sys


def main():
    """Main entry point for the autosarm CLI."""
    parser = argparse.ArgumentParser(
        prog="autosarm",
        description="AutoSARM - Automatic Structure-Activity Relationship Matrix Generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create SAR matrix
  autosarm sarm --csvFile compounds.csv --column IC50_uM --save_folder SAR_Results
  
  # Create SAR tree
  autosarm tree --fragment_core "c1ccccc1" --workFolder SAR_Results --rootTitle Example
  
  # Print SAR table
  autosarm print --core "c1ccc(*)*cc1" --actFile activity.csv

For more information, see: https://github.com/yourorg/autosarm
"""
    )
    
    parser.add_argument(
        "--version", "-v",
        action="version",
        version="%(prog)s 1.0.0"
    )
    
    subparsers = parser.add_subparsers(
        dest="command",
        title="commands",
        description="Available commands"
    )
    
    # SARM subcommand
    sarm_parser = subparsers.add_parser(
        "sarm",
        help="Create SAR matrix from compound library",
        description="Fragment molecules and create SAR matrices"
    )
    _add_sarm_arguments(sarm_parser)
    
    # Tree subcommand
    tree_parser = subparsers.add_parser(
        "tree",
        help="Create SAR tree visualization",
        description="Generate hierarchical SAR tree from existing SAR results"
    )
    _add_tree_arguments(tree_parser)
    
    # Print subcommand
    print_parser = subparsers.add_parser(
        "print",
        help="Print SAR table for a specific core",
        description="Generate and print SAR table for a given core structure"
    )
    _add_print_arguments(print_parser)
    
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        sys.exit(0)
    
    if args.command == "sarm":
        from autosarm.cli.create_sarm import run_sarm
        run_sarm(args)
    elif args.command == "tree":
        from autosarm.cli.create_tree import run_tree
        run_tree(args)
    elif args.command == "print":
        from autosarm.cli.print_sar import run_print
        run_print(args)


def _add_sarm_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for SARM creation command."""
    parser.add_argument(
        "--csvFile",
        required=True,
        help="Path to CSV file with SMILES and activity data"
    )
    parser.add_argument(
        "--type",
        default="smiles",
        choices=["smiles", "scaffold"],
        help="Type of molecular representation (default: smiles)"
    )
    parser.add_argument(
        "--column",
        nargs="+",
        required=True,
        help="Activity column name(s) to include in analysis"
    )
    parser.add_argument(
        "--log",
        type=int,
        default=0,
        help="Apply log transform to activity values (1=yes, 0=no)"
    )
    parser.add_argument(
        "--minimumSite1",
        type=float,
        default=3,
        help="Minimum fragment count for site 1 (default: 3)"
    )
    parser.add_argument(
        "--minimumSite2",
        type=float,
        default=3,
        help="Minimum fragment count for site 2 (default: 3)"
    )
    parser.add_argument(
        "--n_jobs",
        type=int,
        default=8,
        help="Number of parallel jobs (default: 8)"
    )
    parser.add_argument(
        "--save_folder",
        default="SAR_Results",
        help="Output directory for results (default: SAR_Results)"
    )
    parser.add_argument(
        "--csv2excel",
        type=int,
        default=0,
        help="Generate Excel files with molecular images (1=yes, 0=no)"
    )


def _add_tree_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for tree creation command."""
    parser.add_argument(
        "--fragment_core",
        required=True,
        help="Core fragment SMILES for the tree root"
    )
    parser.add_argument(
        "--rootTitle",
        required=True,
        help="Title for the tree (used in output filenames)"
    )
    parser.add_argument(
        "--workFolder",
        required=True,
        help="Working folder containing SAR results"
    )
    parser.add_argument(
        "--inputFile",
        default="input.csv",
        help="Input CSV filename (default: input.csv)"
    )
    parser.add_argument(
        "--maxLevel",
        type=int,
        default=5,
        help="Maximum tree depth (default: 5)"
    )
    parser.add_argument(
        "--treeContent",
        type=str,
        default="['double-cut']",
        help="Tree content types (default: \"['double-cut']\")"
    )
    parser.add_argument(
        "--highlightDict",
        type=str,
        default="[]",
        help="Highlighting configuration as Python dict string"
    )


def _add_print_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for print command."""
    parser.add_argument(
        "--core",
        required=True,
        help="Core structure SMILES with dummy atoms (*)"
    )
    parser.add_argument(
        "--actFile",
        required=True,
        help="Activity data CSV file"
    )
    parser.add_argument(
        "--smiCol",
        default="smiles",
        help="SMILES column name (default: smiles)"
    )
    parser.add_argument(
        "--actCol",
        nargs="+",
        default=["IC50"],
        help="Activity column name(s) (default: IC50)"
    )
    parser.add_argument(
        "--output",
        default="sar_table.csv",
        help="Output filename (default: sar_table.csv)"
    )


if __name__ == "__main__":
    main()
