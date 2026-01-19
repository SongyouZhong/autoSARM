"""
AutoSARM - Automatic Structure-Activity Relationship Matrix Generator

A comprehensive toolkit for drug discovery research that enables automated
SAR analysis and molecular fragmentation.
"""

__version__ = "1.0.0"
__author__ = "Jie Zhang"
__email__ = "jie.zhang@example.com"

from autosarm.core.sarm import create_sarm_matrix, fragmentize
from autosarm.core.tree import create_sar_tree
from autosarm.utils.mol_utils import (
    get_mol,
    canonic_smiles,
    compute_fingerprint,
    compute_similarity,
)

__all__ = [
    "__version__",
    "create_sarm_matrix",
    "fragmentize",
    "create_sar_tree",
    "get_mol",
    "canonic_smiles",
    "compute_fingerprint",
    "compute_similarity",
]
