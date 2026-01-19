"""
CLI module for creating SAR trees.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

from autosarm.core.tree import create_sar_tree


def run_tree(args) -> None:
    """
    Run SAR tree creation.
    
    Args:
        args: Parsed command line arguments
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    
    logging.info("Starting SAR Tree Generation")
    
    # Parse tree content
    try:
        tree_content = eval(args.treeContent)
    except Exception:
        tree_content = ['double-cut']
    
    # Parse highlight dict
    try:
        highlight_dict = eval(args.highlightDict)
    except Exception:
        highlight_dict = []
    
    create_sar_tree(
        fragment_core=args.fragment_core,
        root_title=args.rootTitle,
        work_folder=args.workFolder,
        input_file=args.inputFile,
        tree_content=tree_content,
        highlight_dict=highlight_dict,
        max_level=args.maxLevel
    )
    
    logging.info("SAR Tree Generation Complete")


def get_parser() -> argparse.ArgumentParser:
    """Create argument parser for standalone usage."""
    parser = argparse.ArgumentParser(
        description="Create SAR tree visualization"
    )
    parser.add_argument("--fragment_core", required=True)
    parser.add_argument("--rootTitle", required=True)
    parser.add_argument("--workFolder", required=True)
    parser.add_argument("--inputFile", default="input.csv")
    parser.add_argument("--maxLevel", type=int, default=5)
    parser.add_argument("--treeContent", type=str, default="['double-cut']")
    parser.add_argument("--highlightDict", type=str, default="[]")
    return parser


def main():
    """Main entry point for standalone script."""
    parser = get_parser()
    args = parser.parse_args()
    run_tree(args)


if __name__ == "__main__":
    main()
