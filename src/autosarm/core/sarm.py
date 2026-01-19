"""
Core SARM (Structure-Activity Relationship Matrix) generation module.

This module provides functions for fragmentizing molecules and creating
SAR matrices from compound libraries.
"""

from __future__ import annotations

import copy
import logging
import re
from functools import partial
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs

from autosarm.utils.mol_utils import (
    get_mol,
    safe_mol_from_smiles,
    compute_fingerprint,
    compute_similarity,
    mapper,
)
from autosarm.utils.data_utils import (
    get_number_from_string,
    float_row,
    load_csv_or_df,
    csv_to_excel,
)

logger = logging.getLogger(__name__)

# Try to import pandarallel for parallel processing
try:
    from pandarallel import pandarallel
    HAS_PANDARALLEL = True
except ImportError:
    HAS_PANDARALLEL = False
    logger.warning("pandarallel not available, using sequential processing")


def frags_count(file_or_df: Union[str, pd.DataFrame], minimum_count: int = 5) -> pd.DataFrame:
    """
    Count fragment occurrences and filter by minimum count.
    
    Args:
        file_or_df: CSV file path or DataFrame with 'Key' column
        minimum_count: Minimum occurrence threshold
        
    Returns:
        DataFrame with fragment counts
    """
    df = load_csv_or_df(file_or_df)
    counts = df['Key'].value_counts()
    df_count = pd.DataFrame(counts)
    df_count = df_count[df_count['count'] >= minimum_count]
    df_count['Count'] = df_count['count']
    df_count['Key'] = df_count.index
    df_count = df_count[['Key', 'Count']]
    return df_count


def attach_left_right(frag: str, smi: str) -> Optional[int]:
    """
    Detect if fragment attaches to left (0) or right (1) of molecule.
    
    Args:
        frag: Fragment SMARTS with two dummy atoms
        smi: SMILES of the molecule
        
    Returns:
        0 for left attachment, 1 for right, None on error
    """
    frag_mol = Chem.MolFromSmarts(frag)
    if frag_mol is None:
        return None
        
    match_dummy = frag_mol.GetSubstructMatches(Chem.MolFromSmarts('[#0]'))
    if len(match_dummy) != 2:
        logger.debug(f"Fragment '{frag}' does not have exactly 2 dummy atoms")
        return None
    
    left_dummy = match_dummy[0][0]
    right_dummy = match_dummy[1][0]
    
    mol = safe_mol_from_smiles(smi)
    if mol is None:
        return None
    
    match = mol.GetSubstructMatches(frag_mol)
    if not match:
        return None
        
    atoms = mol.GetAtoms()
    left_atm_idx = match[0][left_dummy]
    right_atm_idx = match[0][right_dummy]
    
    if len(atoms[left_atm_idx].GetNeighbors()) > 1:
        return 0
    if len(atoms[right_atm_idx].GetNeighbors()) > 1:
        return 1
    
    return None


def has_match(smi: str, frag: str = '') -> bool:
    """
    Check if molecule contains a substructure.
    
    Args:
        smi: SMILES string
        frag: SMARTS pattern to match
        
    Returns:
        True if match found
    """
    frag_mol = Chem.MolFromSmarts(frag)
    mol = safe_mol_from_smiles(smi)
    
    if frag_mol is None or mol is None:
        return False
    
    return len(mol.GetSubstructMatches(frag_mol)) > 0


def frag_mol_near_ring(
    smi: str, 
    pos_args: Dict[str, bool] = None
) -> Optional[List[List[str]]]:
    """
    Fragment molecule by cutting bonds near rings.
    
    Args:
        smi: SMILES string to fragment
        pos_args: Dictionary with fragmentation options:
            - 'RR': Cut ring-ring single bonds
            - 'nRnR': Cut all non-ring single bonds
            
    Returns:
        List of [core, fragment, original_smiles] triplets
    """
    if pos_args is None:
        pos_args = {'RR': True, 'nRnR': False}
    
    try:
        mol = get_mol(smi)
        if mol is None:
            return None
        
        if pos_args.get("nRnR", False):
            # All single bonds except those in ring
            bs = []
            for bond in mol.GetBonds():
                idx = bond.GetIdx()
                if bond.IsInRing():
                    continue
                if bond.GetBondType() == Chem.BondType.SINGLE:
                    bs.append(idx)
        else:
            # Bonds connecting non-ring to ring atoms
            bis = mol.GetSubstructMatches(Chem.MolFromSmarts('[!R][R]'))
            bs = [mol.GetBondBetweenAtoms(x, y).GetIdx() for x, y in bis]
            
            if pos_args.get('RR', True):
                # Single bonds between two rings
                bis = mol.GetSubstructMatches(Chem.MolFromSmarts('[R]!@;-[R]'))
                bs.extend([mol.GetBondBetweenAtoms(x, y).GetIdx() for x, y in bis])
        
        frag_pairs = []
        for ibond in bs:
            nm = Chem.FragmentOnBonds(mol, [ibond], dummyLabels=[(0, 0)])
            frags = Chem.GetMolFrags(nm, asMols=True)
            
            # Skip if dummy atoms are adjacent
            dummy_pair = False
            for frag in frags:
                if frag.GetSubstructMatches(Chem.MolFromSmarts('[#0][#0]')):
                    dummy_pair = True
                    break
            
            if dummy_pair:
                continue
            
            smi0 = Chem.MolToSmiles(frags[0])
            smi1 = Chem.MolToSmiles(frags[1])
            
            # Put larger fragment first (as core)
            if (frags[0].GetNumAtoms() < frags[1].GetNumAtoms() and 
                len(frags[0].GetSubstructMatches(Chem.MolFromSmarts('[#0]'))) < 2):
                frag_pairs.append([smi1, smi0, smi])
            else:
                frag_pairs.append([smi0, smi1, smi])
        
        return frag_pairs
        
    except Exception as e:
        logger.debug(f"Error fragmenting molecule: {e}")
        return None


def fragmentize(
    act_cpds: Union[List[str], pd.Series],
    n_jobs: int = 20,
    drop_duplicate: bool = True,
    pos_args: Dict[str, bool] = None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Perform two-round molecular fragmentation.
    
    First round: Fragment all molecules
    Second round: Further fragment the core structures
    
    Args:
        act_cpds: List of SMILES strings
        n_jobs: Number of parallel jobs
        drop_duplicate: Remove duplicate fragments
        pos_args: Fragmentation options
        
    Returns:
        Tuple of (round1_df, round2_df) with fragmentation results
    """
    if pos_args is None:
        pos_args = {'RR': True, 'nRnR': False}
    
    frag_func = partial(frag_mol_near_ring, pos_args=pos_args)
    
    # First round fragmentation
    frag_pairs = mapper(n_jobs)(frag_func, act_cpds)
    frag_pairs = [fp for fp in frag_pairs if fp is not None]
    frag_pairs = [f for fp in frag_pairs for f in fp]
    
    df_round1 = pd.DataFrame(frag_pairs, columns=['Key', 'Value', 'OrgSmi'])
    
    # Second round fragmentation
    if drop_duplicate:
        df_round1_unique = df_round1.drop_duplicates('Key')
    else:
        df_round1_unique = df_round1
    
    frag_pairs2 = mapper(n_jobs)(frag_func, df_round1_unique['Key'])
    frag_pairs2 = [fp for fp in frag_pairs2 if fp is not None]
    frag_pairs2 = [f for fp in frag_pairs2 for f in fp]
    
    df_round2 = pd.DataFrame(frag_pairs2, columns=['Key', 'Value', 'OrgSmi'])
    
    return df_round1, df_round2


def sort_mol_sim(smi_list: Union[List[str], pd.Series]) -> pd.DataFrame:
    """
    Sort SMILES by molecular similarity.
    
    Args:
        smi_list: List of SMILES strings
        
    Returns:
        DataFrame with sorted SMILES
    """
    df_sorted = pd.DataFrame(columns=['Smi', 'Mol', 'Fp'])
    
    for smi in smi_list:
        mol = get_mol(smi)
        fp = compute_fingerprint(mol)
        
        if len(df_sorted) == 0:
            df_sorted.loc[0] = [smi, mol, fp]
            continue
        
        sim_array = np.array([
            DataStructs.TanimotoSimilarity(fp, existing_fp)
            for existing_fp in df_sorted['Fp']
        ])
        max_idx = np.argmax(sim_array)
        
        df_sorted.loc[max_idx + 0.5] = [smi, mol, fp]
        df_sorted = df_sorted.sort_index().reset_index(drop=True)
    
    return pd.DataFrame(df_sorted['Smi'])


def sort_mol_sim_df(df: pd.DataFrame, col: str = 'OrgSmi') -> pd.DataFrame:
    """
    Sort DataFrame by molecular similarity.
    
    Args:
        df: Input DataFrame
        col: Column containing SMILES
        
    Returns:
        Reordered DataFrame
    """
    df = df.copy()
    df['index'] = df[col]
    df = df.set_index('index')
    
    df_sorted = pd.DataFrame(columns=['Key2', 'Value2', 'Smi', 'Mol', 'Fp'])
    
    for idx, row in df.iterrows():
        smi = row[col]
        key2 = row.get('Key', '')
        value2 = row.get('Value', '')
        mol = get_mol(smi)
        fp = compute_fingerprint(mol)
        
        if len(df_sorted) == 0:
            df_sorted.loc[0] = [key2, value2, smi, mol, fp]
            continue
        
        sim_array = np.array([
            DataStructs.TanimotoSimilarity(fp, existing_fp)
            for existing_fp in df_sorted['Fp']
        ])
        max_idx = np.argmax(sim_array)
        
        df_sorted.loc[max_idx + 0.5] = [key2, value2, smi, mol, fp]
        df_sorted = df_sorted.sort_index().reset_index(drop=True)
    
    return pd.DataFrame(df_sorted[['Key2', 'Value2', 'Smi']])


def key_filter(
    key: str,
    ring: bool = True,
    num_atoms: int = 6,
    aromatic: bool = False
) -> bool:
    """
    Filter fragment keys based on structural criteria.
    
    Args:
        key: Fragment SMILES
        ring: Require ring structure
        num_atoms: Minimum number of atoms
        aromatic: Require aromatic ring
        
    Returns:
        True if fragment passes filter
    """
    mol = get_mol(key)
    if mol is None:
        return False
    
    ring_info = mol.GetRingInfo()
    
    if mol.GetNumAtoms() < num_atoms:
        return False
    
    if ring:
        if len(ring_info.AtomRings()) < 1:
            return False
        
        if aromatic:
            aromatic_found = False
            for ring_atoms in ring_info.AtomRings():
                for atom_idx in ring_atoms:
                    if mol.GetAtomWithIdx(atom_idx).GetIsAromatic():
                        aromatic_found = True
                        break
                if aromatic_found:
                    break
            
            if not aromatic_found:
                return False
    
    return True


def replace_nH(ismarts: str, nH: bool = True, DH: bool = True) -> str:
    """
    Normalize nitrogen ionization states in SMARTS.
    
    Args:
        ismarts: Input SMARTS string
        nH: Replace [nH] with n
        DH: Remove deuterium
        
    Returns:
        Normalized SMARTS
    """
    ismarts = ismarts.replace('-', '')
    ismarts = ismarts.replace('\\', '')
    ismarts = ismarts.replace('/', '')
    
    if DH:
        pattern = re.compile(r'\(\[2H\]\)|\[2H\]|\(\)')
        ismarts = re.sub(pattern, '', ismarts)
    
    if nH:
        ismarts = ismarts.replace('[nH]', 'n')
    
    return ismarts


def match_frag(smi: str, ismarts: str) -> int:
    """
    Check if molecule matches fragment pattern.
    
    Args:
        smi: SMILES string
        ismarts: SMARTS pattern
        
    Returns:
        1 if match, 0 otherwise
    """
    mol = safe_mol_from_smiles(smi)
    if mol is None:
        return 0
    
    ismarts = replace_nH(ismarts)
    smarts_mol = Chem.MolFromSmarts(ismarts)
    
    if smarts_mol is None:
        return 0
    
    matched = mol.GetSubstructMatches(smarts_mol)
    return 1 if matched else 0


def set_isotope(mol: Chem.Mol, isotope: int) -> Chem.Mol:
    """
    Set isotope on first dummy atom.
    
    Args:
        mol: RDKit Mol object
        isotope: Isotope number
        
    Returns:
        Modified molecule
    """
    match_list = mol.GetSubstructMatches(Chem.MolFromSmarts('[#0]'))
    if match_list:
        atoms = mol.GetAtoms()
        dummy_idx = match_list[0][0]
        atoms[dummy_idx].SetIsotope(isotope)
    return mol


def get_dummy_neighbor(atom) -> int:
    """Get neighbor index of dummy atom."""
    neighbor = atom.GetNeighbors()[0]
    return int(neighbor.GetIdx())


def connect_R1(
    R: str,
    core: str = '',
    left_or_right: int = 0,
    return_type: str = 'smiles'
) -> Optional[Union[str, Chem.Mol]]:
    """
    Connect an R-group to a core structure.
    
    Args:
        R: R-group SMILES with dummy atom
        core: Core SMILES with two dummy atoms
        left_or_right: Which dummy to connect (0=left, 1=right)
        return_type: 'smiles' or 'mol'
        
    Returns:
        Connected molecule as SMILES or Mol object
    """
    core_mol = get_mol(core)
    if core_mol is None:
        return None
    
    core_mol = copy.deepcopy(core_mol)
    
    match_dummy = core_mol.GetSubstructMatches(Chem.MolFromSmarts('[#0]'))
    if len(match_dummy) != 2:
        logger.debug(f"Core '{core}' does not have exactly 2 dummy atoms")
        return None
    
    left_dummy = match_dummy[0][0]
    right_dummy = match_dummy[1][0]
    
    core_atoms = core_mol.GetAtoms()
    core_atoms[left_dummy].SetIsotope(0)
    core_atoms[right_dummy].SetIsotope(1)
    
    r_mol = safe_mol_from_smiles(R)
    if r_mol is None:
        return None
    
    r_mol = set_isotope(r_mol, 2)
    
    combo = Chem.CombineMols(core_mol, r_mol)
    match = combo.GetSubstructMatches(Chem.MolFromSmarts('[#0]'))
    combo_atoms = combo.GetAtoms()
    
    dummy_pair = []
    dummy_neighbor = []
    
    for m in match:
        atm_idx = m[0]
        isotope = combo_atoms[atm_idx].GetIsotope()
        
        if isotope in [2, left_or_right]:
            dummy_pair.append(atm_idx)
            dummy_neighbor.append(get_dummy_neighbor(combo_atoms[atm_idx]))
        else:
            combo_atoms[atm_idx].SetAtomicNum(50)  # Protect with Sn
    
    if len(dummy_neighbor) < 2:
        return None
    
    ed_combo = Chem.EditableMol(combo)
    ed_combo.AddBond(dummy_neighbor[0], dummy_neighbor[1], order=Chem.rdchem.BondType.SINGLE)
    
    dummy_pair.sort(reverse=True)
    for dummy_idx in dummy_pair:
        ed_combo.RemoveAtom(dummy_idx)
    
    combo = ed_combo.GetMol()
    
    # Replace remaining dummies with H
    products = Chem.ReplaceSubstructs(
        combo, Chem.MolFromSmarts('[#0]'), Chem.MolFromSmarts('[#1]'), replaceAll=True
    )
    combo = products[0]
    
    # Replace Sn back to dummy
    products = Chem.ReplaceSubstructs(
        combo, Chem.MolFromSmarts('[#50]'), Chem.MolFromSmarts('[#0]'), replaceAll=True
    )
    combo = products[0]
    
    combo_smi = Chem.MolToSmiles(combo)
    combo = safe_mol_from_smiles(combo_smi)
    
    if combo is None:
        return None
    
    combo = Chem.RemoveHs(combo)
    
    if return_type == 'mol':
        return combo
    return Chem.MolToSmiles(combo)


def get_activity_info(
    df_act: pd.DataFrame,
    frag: str,
    actCols: List[str],
    smilesCol: str = 'smiles'
) -> str:
    """
    Get activity statistics for a fragment.
    
    Args:
        df_act: DataFrame with activity data
        frag: Fragment SMILES
        actCols: Activity column names
        smilesCol: SMILES column name
        
    Returns:
        Formatted activity string
    """
    df_act = df_act.copy()
    df_act['similarity'] = df_act.apply(
        lambda x: compute_similarity(x[smilesCol], frag),
        axis=1
    )
    
    df_match = df_act[df_act['similarity'] > 0.99]
    df_match = df_match[actCols]
    
    stats = {'mean': {}, 'std': {}, 'median': {}, 'min': {}, 'max': {}}
    
    for col in actCols:
        df_col = df_match.dropna(subset=[col])
        
        if len(df_col) == 0:
            continue
        elif len(df_col) == 1:
            val = float(df_col[col].iloc[0])
            stats['mean'][col] = val
            stats['std'][col] = val
            stats['median'][col] = val
            stats['min'][col] = val
            stats['max'][col] = val
        else:
            stats['mean'][col] = float(df_col[col].mean())
            stats['std'][col] = float(df_col[col].std())
            stats['median'][col] = float(df_col[col].median())
            stats['min'][col] = float(df_col[col].min())
            stats['max'][col] = float(df_col[col].max())
    
    # Format output
    act_info = ''
    for col in actCols:
        if col not in stats['mean']:
            continue
        
        single_value = (len(df_match.dropna(subset=[col])) == 1)
        info = '{0: <15}'.format(f"{col}:")
        
        for metric in ['mean', 'std', 'median', 'min', 'max']:
            if col in stats[metric]:
                info += '{0: <10}'.format(f"{round(stats[metric][col], 1)}")
            if single_value:
                break
        
        info += "|\n"
        act_info += info
    
    return act_info


def create_sarm_matrix(
    df_frags_round1: pd.DataFrame,
    df_frags_round2: pd.DataFrame,
    df_active: pd.DataFrame,
    save_folder: str,
    smi_col: str = "SMILES",
    value_col: List[str] = None,
    minimum_count_site1: int = 5,
    minimum_count_site2: int = 5,
    csv2excel: bool = True,
    cal_table_stats: bool = True,
    n_jobs: int = 50
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Create SAR matrices from fragmentation data.
    
    This is the main function for generating SAR tables including:
    - Left tables (left attachment site variations)
    - Right tables (right attachment site variations)  
    - Combined tables (both sites)
    - Single-cut tables
    
    Args:
        df_frags_round1: First round fragmentation results
        df_frags_round2: Second round fragmentation results
        df_active: Activity data DataFrame
        save_folder: Output directory
        smi_col: Column name for SMILES
        value_col: Activity value column names
        minimum_count_site1: Minimum fragment count for site 1
        minimum_count_site2: Minimum fragment count for site 2
        csv2excel: Generate Excel files with images
        cal_table_stats: Calculate statistics per table
        n_jobs: Number of parallel jobs
        
    Returns:
        Tuple of (left_info_df, right_info_df, combined_info_df)
    """
    if value_col is None:
        value_col = ['Value']
    
    # Initialize pandarallel if available
    if HAS_PANDARALLEL:
        pandarallel.initialize(nb_workers=n_jobs, progress_bar=False)
    
    # Create output directories
    root_path = Path(save_folder).absolute()
    root_path.joinpath('Left_Table').mkdir(exist_ok=True, parents=True)
    root_path.joinpath('Right_Table').mkdir(exist_ok=True, parents=True)
    root_path.joinpath('Combine_Table').mkdir(exist_ok=True, parents=True)
    
    # Prepare activity data
    df_active = df_active.set_index(smi_col, drop=False)
    for col in value_col:
        df_active[col] = df_active[col].apply(get_number_from_string)
    
    # Process round 1 fragments
    logger.info("Processing round 1 fragments")
    frag_round1 = load_csv_or_df(df_frags_round1)
    frag_round1['index'] = frag_round1['Key']
    frag_round1 = frag_round1.set_index('index')
    
    frag_round1_count = frags_count(df_frags_round1, minimum_count=minimum_count_site1)
    frag_round1_count.to_csv(root_path / 'Frag_round1_count.csv', index=None)
    
    if csv2excel:
        csv_to_excel(
            str(root_path / 'Frag_round1_count.csv'),
            img_cols=['Key'],
            save_file=str(root_path / 'Frag_round1_count.xlsx')
        )
    
    # Process round 2 fragments
    logger.info("Processing round 2 fragments")
    frag_round2 = load_csv_or_df(df_frags_round2)
    frag_round2_count = frags_count(df_frags_round2, minimum_count=minimum_count_site2)
    
    # Sort by similarity
    index_list = sort_mol_sim(frag_round2_count['Key'])
    frag_round2_count['index'] = frag_round2_count['Key']
    frag_round2_count = frag_round2_count.set_index('index')
    frag_round2_count = frag_round2_count.loc[index_list['Smi']]
    frag_round2_count = frag_round2_count.reset_index(drop=True)
    frag_round2_count.to_csv(root_path / 'Frag_round2_count.csv', index=None)
    
    if csv2excel:
        csv_to_excel(
            str(root_path / 'Frag_round2_count.csv'),
            img_cols=['Key'],
            save_file=str(root_path / 'Frag_round2_count.xlsx')
        )
    
    # Initialize info DataFrames
    df_table_info_left = pd.DataFrame(columns=['Table', 'Key2', 'Size', 'Items_count'])
    df_table_info_right = pd.DataFrame(columns=['Table', 'Key2', 'Size', 'Items_count'])
    df_table_info_combine = pd.DataFrame(columns=['Table', 'Key2', 'Size', 'Items_count'])
    df_table_info_singlecut = pd.DataFrame(columns=['Table', 'Key2', 'Size', 'Items_count'])
    
    # Create single-cut tables
    logger.info("Creating single-cut tables")
    frag_round1_count = frag_round1_count.reset_index(drop=True)
    
    for idx, row in frag_round1_count.iterrows():
        try:
            logger.debug(f"Processing single-cut row: {idx}")
            round1_key = row['Key']
            round1_key_list = frag_round1[frag_round1['Key'] == round1_key]
            
            df_round1_key = pd.DataFrame(round1_key_list)
            df_round1_key_sort = sort_mol_sim_df(df_round1_key, col='OrgSmi')
            
            df_sarm = df_round1_key_sort.copy()
            df_sarm['index'] = df_sarm['Smi']
            df_sarm = df_sarm.set_index('index')
            df_sarm = df_sarm.reset_index(drop=True)
            df_sarm.loc[-0.5] = list(df_sarm.count(axis=0))
            df_sarm = df_sarm.sort_index()
            df_sarm.insert(1, 'Count', list(df_sarm.count(axis=1) - 3))
            
            df_tmp = pd.DataFrame({
                'Table': [f'Table_{idx}_singleCut.csv'],
                'Key2': [round1_key],
                'Size': [df_sarm.shape],
                'Items_count': [df_sarm.iloc[0].astype(int)[3:].sum()]
            })
            
            df_sarm = df_sarm.drop_duplicates('Value2')
            df_sarm = df_sarm.reset_index(drop=True)
            
            # Add activity information
            for jdx, jrow in df_sarm.iterrows():
                if jrow['Smi'] not in df_active.index:
                    continue
                df_sarm.loc[jdx, 'activity'] = get_activity_info(
                    df_active, jrow['Smi'], actCols=value_col, smilesCol='smiles'
                )
            
            df_sarm.to_csv(root_path / f'Combine_Table/Table_{idx}_singleCut.csv', index=None)
            
            if csv2excel:
                csv_to_excel(
                    str(root_path / f'Combine_Table/Table_{idx}_singleCut.csv'),
                    img_cols=['Value2', 'Smi'],
                    save_file=str(root_path / f'Left_Table/Table_{idx}_left.xlsx')
                )
            
            df_table_info_singlecut = pd.concat([df_table_info_singlecut, df_tmp])
            
        except Exception as e:
            logger.warning(f"Error in single-cut row {idx}: {e}")
            continue
    
    df_table_info_singlecut = df_table_info_singlecut.sort_values('Items_count', ascending=False)
    df_table_info_singlecut = df_table_info_singlecut.reset_index(drop=True)
    
    # Create double-cut tables
    logger.info("Creating double-cut tables")
    
    for idx, row in frag_round2_count.iterrows():
        try:
            logger.debug(f"Processing double-cut row: {idx}")
            key2 = row['Key']
            
            if not key_filter(key2, ring=True, num_atoms=6, aromatic=True):
                continue
            
            df_key2 = frag_round2[frag_round2['Key'] == key2]
            df_key2 = pd.DataFrame(df_key2)
            
            break_bond_pos = [
                attach_left_right(r['Key'], r['OrgSmi']) 
                for _, r in df_key2.iterrows()
            ]
            df_key2['BreakBondPos'] = break_bond_pos
            
            df_sarm_list = []
            
            for attach_id in [0, 1]:
                df_left = df_key2[df_key2['BreakBondPos'] == attach_id]
                
                if len(df_left) == 0:
                    df_sarm_list.append(pd.DataFrame())
                    continue
                
                df_key1_sort = sort_mol_sim_df(df_left, col='OrgSmi')
                df_sarm = df_key1_sort.copy()
                df_sarm['index'] = df_sarm['Smi']
                df_sarm = df_sarm.set_index('index')
                
                # Calculate statistics
                if cal_table_stats and len(df_key1_sort) > 0:
                    def calc_stats(irow):
                        df_key1_sel = frag_round1[frag_round1['Key'] == irow['Smi']]
                        for _, jrow in df_key1_sel.iterrows():
                            return [
                                irow['Smi'],
                                str(jrow['Value']),
                                get_activity_info(df_active, jrow['OrgSmi'], actCols=value_col, smilesCol='smiles')
                            ]
                        return [irow['Smi'], '', '']
                    
                    if HAS_PANDARALLEL:
                        stats_list = df_key1_sort.parallel_apply(calc_stats, axis=1)
                    else:
                        stats_list = df_key1_sort.apply(calc_stats, axis=1)
                    
                    for stat in stats_list:
                        if stat[1]:
                            df_sarm.loc[stat[0], stat[1]] = stat[2]
                
                df_sarm = df_sarm.reset_index(drop=True)
                df_sarm.loc[-0.5] = list(df_sarm.count(axis=0))
                df_sarm = df_sarm.sort_index()
                df_sarm.insert(1, 'Count', list(df_sarm.count(axis=1) - 3))
                
                df_tmp = pd.DataFrame({
                    'Table': [f'Table_{idx}_{"left" if attach_id == 0 else "right"}.csv'],
                    'Key2': [key2],
                    'Size': [df_sarm.shape],
                    'Items_count': [df_sarm.iloc[0].astype(int)[4:].sum() if len(df_sarm) > 0 else 0]
                })
                
                df_sarm = df_sarm.drop_duplicates('Value2')
                
                if attach_id == 0:
                    df_sarm.to_csv(root_path / f'Left_Table/Table_{idx}_left.csv', index=None)
                    if csv2excel:
                        csv_to_excel(
                            str(root_path / f'Left_Table/Table_{idx}_left.csv'),
                            img_cols=['Value2', 'Smi'],
                            save_file=str(root_path / f'Left_Table/Table_{idx}_left.xlsx')
                        )
                    df_table_info_left = pd.concat([df_table_info_left, df_tmp])
                else:
                    df_sarm.to_csv(root_path / f'Right_Table/Table_{idx}_right.csv', index=None)
                    if csv2excel:
                        csv_to_excel(
                            str(root_path / f'Right_Table/Table_{idx}_right.csv'),
                            img_cols=['Value2', 'Smi'],
                            save_file=str(root_path / f'Right_Table/Table_{idx}_right.xlsx')
                        )
                    df_table_info_right = pd.concat([df_table_info_right, df_tmp])
                
                df_sarm_list.append(df_sarm.copy())
            
            # Combine left and right tables
            if len(df_sarm_list) == 2 and len(df_sarm_list[0]) > 0 and len(df_sarm_list[1]) > 0:
                if df_sarm_list[0].size > df_sarm_list[1].size:
                    df_sarm_main = df_sarm_list[0]
                    df_sarm_sub = df_sarm_list[1]
                    connect_site = 0
                else:
                    df_sarm_main = df_sarm_list[1]
                    df_sarm_sub = df_sarm_list[0]
                    connect_site = 1
                
                frag_cols = df_sarm_sub.columns[4:]
                add_frag_cols = [f for f in frag_cols if f not in list(df_sarm_main['Value2'])]
                
                df_tmp_add = pd.DataFrame({'Value2': add_frag_cols})
                df_sarm_main = pd.concat([df_sarm_main, df_tmp_add])
                
                df_sarm_main['index'] = df_sarm_main['Value2']
                df_sarm_main = df_sarm_main.set_index('index')
                
                df_sarm_sub['index'] = df_sarm_sub['Value2']
                df_sarm_sub = df_sarm_sub.set_index('index')
                
                for iidx, sub_row in df_sarm_sub.iterrows():
                    if isinstance(iidx, str):
                        for icol in frag_cols:
                            i_key2 = df_sarm_sub.loc[sub_row['Value2'], 'Key2']
                            df_sarm_main.loc[icol, 'Key2'] = i_key2
                            i_value2 = df_sarm_sub.loc[sub_row['Value2'], icol]
                            df_sarm_main.loc[icol, sub_row['Value2']] = i_value2
                            df_sarm_main.loc[icol, 'Smi'] = connect_R1(
                                R=icol, core=i_key2, left_or_right=connect_site, return_type='smiles'
                            )
                
                df_sarm_main = df_sarm_main.reset_index(drop=True)
                df_sarm_main.loc[0] = list(df_sarm_main.count(axis=0))
                df_sarm_main = df_sarm_main.sort_index()
                
                if 'Count' in df_sarm_main.columns:
                    df_sarm_main = df_sarm_main.drop('Count', axis=1)
                df_sarm_main.insert(1, 'Count', list(df_sarm_main.count(axis=1) - 3))
                
                df_sarm_main.to_csv(root_path / f'Combine_Table/Table_{idx}_combine.csv', index=None)
                
                if csv2excel:
                    csv_to_excel(
                        str(root_path / f'Combine_Table/Table_{idx}_combine.csv'),
                        img_cols=['Value2', 'Smi'],
                        save_file=str(root_path / f'Combine_Table/Table_{idx}_combine.xlsx')
                    )
                
                df_tmp = pd.DataFrame({
                    'Table': [f'Table_{idx}_combine.csv'],
                    'Key2': [key2],
                    'Size': [df_sarm_main.shape],
                    'Items_count': [df_sarm_main.iloc[0].astype(int)[3:].sum() if len(df_sarm_main) > 0 else 0]
                })
                df_table_info_combine = pd.concat([df_table_info_combine, df_tmp])
                
        except Exception as e:
            logger.warning(f"Error in double-cut row {idx}: {e}")
            continue
    
    # Sort and save info tables
    df_table_info_left = df_table_info_left.sort_values('Items_count', ascending=False)
    df_table_info_right = df_table_info_right.sort_values('Items_count', ascending=False)
    df_table_info_combine = df_table_info_combine.sort_values('Items_count', ascending=False)
    df_table_info_combine = df_table_info_combine.reset_index(drop=True)
    
    df_table_info_left.to_csv(root_path / 'Left_Table_info.csv', index=None)
    df_table_info_right.to_csv(root_path / 'Right_Table_info.csv', index=None)
    df_table_info_combine.to_csv(root_path / 'Combine_Table_info.csv', index=None)
    df_table_info_singlecut.to_csv(root_path / 'singleCut_Table_info.csv', index=None)
    
    if csv2excel:
        for fname in ['Left_Table_info', 'Right_Table_info', 'Combine_Table_info', 'singleCut_Table_info']:
            csv_to_excel(
                str(root_path / f'{fname}.csv'),
                img_cols=['Key2'],
                save_file=str(root_path / f'{fname}.xlsx')
            )
    
    logger.info(f"SAR matrix creation complete. Results saved to: {root_path}")
    
    return df_table_info_left, df_table_info_right, df_table_info_combine
