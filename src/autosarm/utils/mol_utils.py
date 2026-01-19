"""
Molecular utilities for AutoSARM.

This module provides core molecular processing functions using RDKit.
"""

from __future__ import annotations

import copy
import logging
import re
from functools import partial
from multiprocessing import Pool
from typing import List, Optional, Union

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, rdFingerprintGenerator

logger = logging.getLogger(__name__)

# Global Morgan fingerprint generator for performance
_morgan_generator: Optional[rdFingerprintGenerator.FingerprintGenerator] = None


def get_morgan_generator(radius: int = 2, fp_size: int = 1024):
    """Get or create a cached MorganGenerator instance for performance."""
    global _morgan_generator
    if _morgan_generator is None:
        _morgan_generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius, fpSize=fp_size
        )
    return _morgan_generator


def is_valid_smiles(smi: Union[str, None]) -> bool:
    """
    Validate if a SMILES string is valid.
    
    Args:
        smi: SMILES string to validate
        
    Returns:
        True if valid, False otherwise
    """
    if not isinstance(smi, str):
        return False
    
    smi_stripped = smi.strip()
    if not smi_stripped or smi_stripped.isdigit():
        return False
    
    try:
        mol = Chem.MolFromSmiles(smi_stripped)
        return mol is not None
    except Exception:
        return False


def get_mol(smiles_or_mol: Union[str, Chem.Mol, None]) -> Optional[Chem.Mol]:
    """
    Load SMILES string or molecule into RDKit Mol object.
    
    Args:
        smiles_or_mol: SMILES string or RDKit Mol object
        
    Returns:
        RDKit Mol object or None if invalid
    """
    if smiles_or_mol is None:
        return None
        
    if isinstance(smiles_or_mol, Chem.Mol):
        return smiles_or_mol
        
    if isinstance(smiles_or_mol, str):
        if len(smiles_or_mol) == 0:
            return None
        try:
            mol = Chem.MolFromSmiles(smiles_or_mol)
            if mol is not None:
                Chem.SanitizeMol(mol)
            return mol
        except (ValueError, RuntimeError):
            return None
    
    return None


def safe_mol_from_smiles(
    smi: Union[str, int, float, None], 
    log_errors: bool = False
) -> Optional[Chem.Mol]:
    """
    Safely create a molecule from SMILES with comprehensive type checking.
    
    Args:
        smi: SMILES string (or numeric value that will be rejected)
        log_errors: Whether to log parsing errors
        
    Returns:
        RDKit Mol object or None if parsing fails
    """
    # Handle numeric types
    if isinstance(smi, (int, float, np.integer, np.floating)):
        if log_errors:
            logger.warning(f"Received numeric value {smi} instead of SMILES string")
        return None
    
    smi_str = str(smi) if smi is not None else ""
    
    if not is_valid_smiles(smi_str):
        if log_errors:
            logger.warning(f"Invalid SMILES: {smi_str}")
        return None
    
    try:
        return Chem.MolFromSmiles(smi_str)
    except Exception as e:
        if log_errors:
            logger.error(f"Error parsing SMILES '{smi_str}': {e}")
        return None


def canonic_smiles(smiles_or_mol: Union[str, Chem.Mol]) -> Optional[str]:
    """
    Convert molecule to canonical SMILES.
    
    Args:
        smiles_or_mol: SMILES string or RDKit Mol object
        
    Returns:
        Canonical SMILES string or None if conversion fails
    """
    try:
        mol = get_mol(smiles_or_mol)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, isomericSmiles=False)
    except Exception as e:
        logger.debug(f"Failed to canonicalize SMILES: {e}")
        return None


def compute_fingerprint(
    mol: Union[str, Chem.Mol], 
    radius: int = 2, 
    n_bits: int = 1024
) -> Optional[DataStructs.ExplicitBitVect]:
    """
    Compute Morgan fingerprint for a molecule.
    
    Args:
        mol: SMILES string or RDKit Mol object
        radius: Fingerprint radius (default: 2)
        n_bits: Number of bits (default: 1024)
        
    Returns:
        Morgan fingerprint or None if computation fails
    """
    mol = get_mol(mol)
    if mol is None:
        return None
    
    try:
        generator = get_morgan_generator(radius=radius, fp_size=n_bits)
        return generator.GetFingerprint(mol)
    except Exception as e:
        logger.warning(f"Error computing fingerprint: {e}")
        return None


# Alias for backward compatibility
compute_FP = compute_fingerprint


def compute_similarity(
    smi: str, 
    smi_list: Union[str, List[str], List[DataStructs.ExplicitBitVect]], 
    mode: str = 'smi-smi'
) -> Union[float, List[float]]:
    """
    Compute Tanimoto similarity between molecules.
    
    Args:
        smi: Query SMILES string
        smi_list: Target SMILES string(s) or fingerprints
        mode: Comparison mode ('smi-smi', 'smi-smis', 'smi-FPs')
        
    Returns:
        Similarity score(s)
    """
    generator = get_morgan_generator(radius=2, fp_size=1024)
    
    if mode == 'smi-smis':
        mol1 = safe_mol_from_smiles(smi)
        if mol1 is None:
            return [0.0] * len(smi_list)
        fp1 = generator.GetFingerprint(mol1)
        
        similarities = []
        for target_smi in smi_list:
            mol2 = safe_mol_from_smiles(target_smi)
            if mol2 is not None:
                fp2 = generator.GetFingerprint(mol2)
                similarities.append(DataStructs.TanimotoSimilarity(fp1, fp2))
            else:
                similarities.append(0.0)
        return similarities
        
    elif mode == 'smi-smi':
        mol1 = safe_mol_from_smiles(smi)
        mol2 = safe_mol_from_smiles(smi_list)
        if mol1 is None or mol2 is None:
            return 0.0
        
        fp1 = generator.GetFingerprint(mol1)
        fp2 = generator.GetFingerprint(mol2)
        return DataStructs.TanimotoSimilarity(fp1, fp2)
        
    elif mode == 'smi-FPs':
        mol1 = safe_mol_from_smiles(smi)
        if mol1 is None:
            return [0.0] * len(smi_list)
        fp1 = generator.GetFingerprint(mol1)
        return [DataStructs.TanimotoSimilarity(fp, fp1) for fp in smi_list]
    
    return 0.0


# Alias for backward compatibility
compute_sim = compute_similarity


def kekulize_smi(smi: str) -> str:
    """
    Convert SMILES to Kekulized form.
    
    Args:
        smi: Input SMILES string
        
    Returns:
        Kekulized SMILES string
    """
    try:
        mol = get_mol(smi)
        if mol is None:
            return smi
        Chem.Kekulize(mol, clearAromaticFlags=True)
        return Chem.MolToSmiles(mol, kekuleSmiles=True)
    except Exception as e:
        logger.debug(f"Kekulization failed: {e}")
        return smi


def heavy_atom_count(smi: str) -> int:
    """
    Count heavy atoms in a molecule.
    
    Args:
        smi: SMILES string
        
    Returns:
        Number of heavy atoms
    """
    heavy_atoms = {'C', 'c', 'N', 'n', 'O', 'o', 'S', 's', 'F', 'Br', 'Cl', 'I', 'B', 'P'}
    return sum(1 for char in smi if char in heavy_atoms)


def remove_dummy(smi: str) -> str:
    """
    Remove dummy atoms (*) from a SMILES string.
    
    Args:
        smi: SMILES string with dummy atoms
        
    Returns:
        SMILES string with dummy atoms replaced by H and removed
    """
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return smi
            
        # Replace [*]-N with H-N
        matched = mol.GetSubstructMatches(Chem.MolFromSmarts("[#7][#0]"))
        atoms = mol.GetAtoms()
        for match in matched:
            atoms[match[1]].SetAtomicNum(1)
        
        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)
        smi_clean = Chem.MolToSmiles(mol)
        
        # Remove remaining dummy atoms
        pattern = re.compile(r'\(\*\)|\*|\[\*\]|\-')
        return re.sub(pattern, '', smi_clean)
    except Exception:
        return smi


def mol_with_atom_index(
    mol: Union[str, Chem.Mol], 
    idx: Optional[int] = None, 
    start: int = 1
) -> Chem.Mol:
    """
    Add atom index labels to molecule for visualization.
    
    Args:
        mol: SMILES string or RDKit Mol object
        idx: Specific index to label all atoms (or None for sequential)
        start: Starting index number (default: 1)
        
    Returns:
        RDKit Mol object with atom notes
    """
    mol = get_mol(mol)
    mol = copy.deepcopy(mol)
    
    for atom in mol.GetAtoms():
        if idx is not None:
            atom.SetProp("atomNote", str(idx))
        else:
            atom.SetProp("atomNote", str(atom.GetIdx() + start))
    
    return mol


def mapper(n_jobs: int):
    """
    Return a map function for parallel processing.
    
    Args:
        n_jobs: Number of parallel jobs (1 for sequential)
        
    Returns:
        Map function
    """
    if n_jobs == 1:
        def _mapper(*args, **kwargs):
            return list(map(*args, **kwargs))
        return _mapper
    
    if isinstance(n_jobs, int):
        pool = Pool(n_jobs)
        
        def _mapper(*args, **kwargs):
            try:
                result = pool.map(*args, **kwargs)
            finally:
                pool.terminate()
            return result
        
        return _mapper
    
    return n_jobs.map


def show_mols(
    smiles_mols: List[Union[str, Chem.Mol]], 
    sub_img_size: tuple = (500, 400), 
    legends: Optional[List[str]] = None,
    mols_per_row: int = 3
):
    """
    Display multiple molecules in a grid.
    
    Args:
        smiles_mols: List of SMILES strings or Mol objects
        sub_img_size: Size of each molecule image
        legends: Optional legends for each molecule
        mols_per_row: Number of molecules per row
        
    Returns:
        Tuple of (PNG image, SVG image)
    """
    legends = legends or []
    mols = [get_mol(smi) for smi in smiles_mols]
    
    mol_cls = []
    legends_cls = []
    
    for i, mol in enumerate(mols):
        if mol is None:
            continue
        mol_cls.append(mol)
        legends_cls.append(legends[i] if i < len(legends) else '')
    
    svg = Draw.MolsToGridImage(
        mol_cls, 
        subImgSize=sub_img_size, 
        molsPerRow=mols_per_row, 
        useSVG=True, 
        legends=legends_cls
    )
    png = Draw.MolsToGridImage(
        mol_cls, 
        subImgSize=sub_img_size, 
        useSVG=False, 
        molsPerRow=mols_per_row, 
        legends=legends_cls
    )
    
    return png, svg
