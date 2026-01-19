"""
SAR Tree generation module.

This module provides functions for creating and visualizing 
Structure-Activity Relationship trees.
"""

from __future__ import annotations

import logging
import os
import re
from functools import partial
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import graphviz
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

from autosarm.utils.mol_utils import (
    get_mol,
    safe_mol_from_smiles,
    canonic_smiles,
    remove_dummy,
    mapper,
)
from autosarm.utils.data_utils import float_row

logger = logging.getLogger(__name__)


class TreeNode:
    """Node in the SAR tree representing a molecular fragment."""
    
    def __init__(self, smiles: str):
        """
        Initialize a tree node.
        
        Args:
            smiles: SMILES string of the fragment
        """
        self.root = True
        self.SMILES = smiles
        self.children: List[TreeNode] = []
        self.parents: List[TreeNode] = []
    
    def add_child(self, node: 'TreeNode') -> None:
        """Add a child node."""
        node.root = False
        self.children.append(node)
    
    def print_tree(
        self, 
        df_table: pd.DataFrame, 
        opfile=None,
        max_levels: int = 100
    ) -> List[List[List[str]]]:
        """
        Print tree structure and return SMILES at each level.
        
        Args:
            df_table: DataFrame with table information
            opfile: Output file handle for text output
            max_levels: Maximum tree depth
            
        Returns:
            Nested list of SMILES strings by level
        """
        tree_list = []
        level = 1
        
        info = fetch_table_info(df_table, self.SMILES)
        print(f"Level {level}:  {self.SMILES} {info}")
        tree_list.append([[self.SMILES]])
        
        if opfile:
            opfile.write(f"Level {level}:  {self.SMILES} {info}\n")
        
        # First level children
        level_smiles = []
        for idx, child in enumerate(self.children):
            child_info = fetch_table_info(df_table, child.SMILES)
            print(f'node_{idx} {child.SMILES} {child_info} |', end=' ')
            if opfile:
                opfile.write(f"node_{idx}  {child.SMILES} {child_info} |")
            level_smiles.append(child.SMILES)
        
        tree_list.append([level_smiles])
        if opfile:
            opfile.write('\n')
        
        # Traverse remaining levels
        next_level_nodes = self.children
        count = 0
        
        while len(next_level_nodes) > 0 and count < max_levels:
            level += 1
            print(f"Level {level}:", end='  ')
            if opfile:
                opfile.write(f"Level {level}:  ")
            
            tmp_next_level = []
            level_smiles = []
            
            for idx, node in enumerate(next_level_nodes):
                print(f'node_{idx}', end=' ')
                if opfile:
                    opfile.write(f"node_{idx}  ")
                
                if node.children:
                    node_smiles = []
                    for child in node.children:
                        child_info = fetch_table_info(df_table, child.SMILES)
                        print(f"{child.SMILES} {child_info}", end=' ')
                        if opfile:
                            opfile.write(f"{child.SMILES} {child_info} ")
                        node_smiles.append(child.SMILES)
                        tmp_next_level.append(child)
                    level_smiles.append(node_smiles)
                else:
                    level_smiles.append([])
                
                print(' |', end=' ')
                if opfile:
                    opfile.write(' |')
            
            tree_list.append(level_smiles)
            print()
            if opfile:
                opfile.write('\n')
            
            next_level_nodes = tmp_next_level
            count += 1
        
        return tree_list


def fetch_table_info(df_table: pd.DataFrame, smi: str) -> str:
    """
    Get table information for a SMILES.
    
    Args:
        df_table: DataFrame with table info indexed by SMILES
        smi: SMILES string to look up
        
    Returns:
        Formatted info string
    """
    try:
        if smi not in df_table.index:
            return ""
        
        info = df_table.loc[smi, ['Size', 'Items_count']]
        
        if isinstance(info['Items_count'], list):
            return f"{info['Size'].tolist()}{info['Items_count'].tolist()}"
        return f"{info['Size']}{info['Items_count']}"
    except Exception:
        return ""


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


def get_atom_ring_info(smi: str) -> np.ndarray:
    """
    Get ring membership info for each atom.
    
    Args:
        smi: SMILES string
        
    Returns:
        Array with ring counts per atom
    """
    mol = get_mol(smi)
    if mol is None:
        return np.array([])
    
    ring_info = mol.GetRingInfo()
    atom_ring_info = np.zeros(mol.GetNumAtoms())
    
    for ring in ring_info.AtomRings():
        atom_ring_info[list(ring)] += 1
    
    return atom_ring_info


def match_core(smi: str, ismarts: str) -> bool:
    """
    Check if molecule matches core structure with exact ring matching.
    
    Args:
        smi: SMILES string
        ismarts: SMARTS pattern
        
    Returns:
        True if exact ring match
    """
    mol = get_mol(smi)
    if mol is None:
        return False
    
    # Remove dummy atoms for matching
    if '*' in ismarts:
        smarts_no_dummy = remove_dummy(ismarts)
    else:
        smarts_no_dummy = ismarts
    
    smarts_mol = Chem.MolFromSmarts(replace_nH(smarts_no_dummy))
    if smarts_mol is None:
        return False
    
    matched = mol.GetSubstructMatches(smarts_mol)
    if not matched:
        return False
    
    smarts_ring_info = get_atom_ring_info(smarts_no_dummy)
    mol_ring_info = get_atom_ring_info(mol)
    
    for match in matched:
        submol_ring_info = mol_ring_info[list(match)]
        diff = np.abs(smarts_ring_info - submol_ring_info).sum()
        if diff == 0:
            return True
    
    return False


def is_parent(parent: str, child: str, ring: bool = False) -> bool:
    """
    Check if parent is a substructure of child.
    
    Args:
        parent: Parent SMARTS
        child: Child SMILES
        ring: Require exact ring matching
        
    Returns:
        True if parent is substructure of child
    """
    mol_parent = Chem.MolFromSmarts(parent)
    mol_child = Chem.MolFromSmiles(child, sanitize=False)
    
    if mol_parent is None or mol_child is None:
        return False
    
    try:
        mol_child.UpdatePropertyCache()
    except Exception:
        return False
    
    matched = mol_child.GetSubstructMatches(mol_parent)
    
    if not matched:
        return False
    
    if ring:
        try:
            return match_core(child, parent)
        except Exception as e:
            logger.debug(f"Ring match error: child={child}, parent={parent}, error={e}")
            return False
    
    return True


def is_child(child: str, parent: str, ring: bool = False) -> bool:
    """Check if child is a superstructure of parent."""
    return is_parent(parent, child, ring)


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


def remove_ionization_form(smi: Union[str, List[str]]) -> Union[str, List[str]]:
    """
    Remove ionization markers from SMILES.
    
    Args:
        smi: Single SMILES or list of SMILES
        
    Returns:
        Cleaned SMILES
    """
    if isinstance(smi, list):
        return [remove_ionization_form(s) for s in smi]
    
    smi = smi.replace('[nH]', 'n')
    smi = smi.replace('[NH]', 'N')
    return smi


def find_children(smi_list: List[str]) -> Dict[str, List[str]]:
    """
    Find all children (superstructures) for each SMILES.
    
    Args:
        smi_list: List of SMILES strings
        
    Returns:
        Dictionary mapping each SMILES to its children
    """
    children_dict = {}
    
    for smi in smi_list:
        children = []
        for other_smi in smi_list:
            if other_smi != smi and is_child(other_smi, smi, ring=True):
                children.append(other_smi)
        children_dict[smi] = children
    
    return children_dict


def find_children_single(smi: str, smi_list: List[str]) -> List[str]:
    """Find children for a single SMILES."""
    children = []
    for other_smi in smi_list:
        if other_smi != smi and is_child(other_smi, smi, ring=True):
            children.append(other_smi)
    return children


def find_parents(smi_list: List[str]) -> Dict[str, List[str]]:
    """
    Find all parents (substructures) for each SMILES.
    
    Args:
        smi_list: List of SMILES strings
        
    Returns:
        Dictionary mapping each SMILES to its parents
    """
    parent_dict = {}
    
    for smi in smi_list:
        parents = []
        for other_smi in smi_list:
            if other_smi != smi and is_parent(other_smi, smi, ring=True):
                parents.append(other_smi)
        parent_dict[smi] = parents
    
    return parent_dict


def find_parents_single(smi: str, smi_list: List[str]) -> List[str]:
    """Find parents for a single SMILES."""
    parents = []
    for other_smi in smi_list:
        if other_smi != smi and is_parent(other_smi, smi, ring=True):
            parents.append(other_smi)
    return parents


def if_root(smi: str, smi_list: List[str]) -> bool:
    """Check if SMILES has no parents in the list."""
    return len(find_parents_single(smi, smi_list)) == 0


def find_root(smi_list: List[str]) -> List[str]:
    """Find all root SMILES (those with no parents)."""
    return [smi for smi in smi_list if if_root(smi, smi_list)]


def real_son_node(
    child_smi: str, 
    children_smis: List[str], 
    parent_dict: Dict[str, List[str]]
) -> bool:
    """
    Check if child is a direct child (not grandchild).
    
    A real child should not have any parents in the children list.
    
    Args:
        child_smi: Potential child SMILES
        children_smis: List of all children
        parent_dict: Parent relationships
        
    Returns:
        True if direct child
    """
    if child_smi not in parent_dict:
        return True
    
    child_parents = parent_dict[child_smi]
    
    for parent in child_parents:
        if parent in children_smis:
            return False
    
    return True


def show_cpd(
    img_cpd: List,
    df_table: pd.DataFrame,
    df_act: pd.DataFrame,
    actCols: List[str],
    highlightDictList: List[Dict]
) -> Optional[List]:
    """
    Generate compound image with activity highlighting.
    
    Args:
        img_cpd: [image_path, smiles] pair
        df_table: Table info DataFrame
        df_act: Activity DataFrame
        actCols: Activity column names
        highlightDictList: Highlighting rules
        
    Returns:
        [image_path, smiles] or None
    """
    try:
        img_path = img_cpd[0]
        cpd_smi = img_cpd[1]
        
        # Get molecule
        mol = safe_mol_from_smiles(cpd_smi)
        if mol is None:
            mol = Chem.MolFromSmarts(cpd_smi)
        
        if mol is None:
            return None
        
        # Get info text
        info = fetch_table_info(df_table, cpd_smi)
        
        # Draw molecule
        drawer = rdMolDraw2D.MolDraw2DCairo(300, 250)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        
        # Save image
        with open(img_path, 'wb') as f:
            f.write(drawer.GetDrawingText())
        
        return [img_path, cpd_smi]
        
    except Exception as e:
        logger.debug(f"Error showing compound: {e}")
        return None


def create_sar_tree(
    fragment_core: str,
    root_title: str,
    work_folder: str,
    input_file: str = 'input.csv',
    tree_content: List[str] = None,
    highlight_dict: List[Dict] = None,
    max_level: int = 5
) -> None:
    """
    Create SAR tree visualization.
    
    Args:
        fragment_core: Core fragment SMILES
        root_title: Title for the tree
        work_folder: Working directory with SAR results
        input_file: Input CSV filename
        tree_content: Content types to include ('double-cut', 'single-cut', 'whole-compound')
        highlight_dict: Highlighting rules for activity
        max_level: Maximum tree depth
    """
    if tree_content is None:
        tree_content = ['double-cut']
    if highlight_dict is None:
        highlight_dict = []
    
    work_folder = os.path.abspath(work_folder)
    
    # Get activity columns from highlight dict
    act_cols = [item['col'] for item in highlight_dict]
    
    # Create tree directory
    tree_path = Path(f"{work_folder}/Trees/FragTree_{root_title}").absolute()
    tree_path.mkdir(exist_ok=True, parents=True)
    
    roots = [fragment_core]
    df_list = []
    
    # Load data based on tree content
    if 'double-cut' in tree_content:
        df_tmp = pd.read_csv(f"{work_folder}/Combine_Table_info.csv")
        if "Key2" in df_tmp.columns:
            df_tmp["SMILES"] = df_tmp["Key2"]
        elif "SMILES" not in df_tmp.columns:
            raise ValueError("Combine_Table_info.csv must have 'Key2' or 'SMILES' column")
        df_list.append(df_tmp)
    
    if 'single-cut' in tree_content:
        df_tmp = pd.read_csv(f"{work_folder}/singleCut_Table_info.csv")
        if "Key2" in df_tmp.columns:
            df_tmp["SMILES"] = df_tmp["Key2"]
        elif "SMILES" not in df_tmp.columns:
            raise ValueError("singleCut_Table_info.csv must have 'Key2' or 'SMILES' column")
        df_list.append(df_tmp)
    
    if 'whole-compound' in tree_content:
        df_tmp = pd.read_csv(f"{work_folder}/{input_file}")
        df_tmp["SMILES"] = df_tmp["smiles"].apply(canonic_smiles)
        df_tmp["Items_count"] = 1
        df_list.append(df_tmp)
    
    df_table = pd.concat(df_list)
    
    # Load activity data
    df_act = pd.read_csv(f"{work_folder}/{input_file}")
    df_act = float_row(df_act, cols=act_cols, dropna=False)
    
    if 'SMILES' not in df_act.columns and 'smiles' in df_act.columns:
        df_act['SMILES'] = df_act['smiles']
    df_act['Cano_SMILES'] = df_act['SMILES'].apply(canonic_smiles)
    df_act = df_act.dropna(subset=['Cano_SMILES'])
    
    # Set up table index
    df_table['Matched'] = df_table.apply(
        lambda x: match_frag(x['SMILES'], ismarts=fragment_core), 
        axis=1
    )
    df_table = df_table[df_table['Matched'] == 1]
    
    smiles_list = df_table["SMILES"].tolist()
    df_table['index'] = df_table["SMILES"]
    df_table.set_index('SMILES', inplace=True, drop=False)
    df_table['Items_count'] = df_table['Items_count'].astype(int)
    
    # Build parent/child relationships
    logger.info("Finding children relationships")
    children_dict = find_children(smiles_list)
    
    logger.info("Finding parent relationships")
    parent_dict = find_parents(smiles_list)
    
    # Build trees
    logger.info("Building tree structure")
    tree_list = []
    tree_smile_list = []
    
    with open(f'{tree_path}/Combine_Table_info_Tree_{root_title}.txt', 'w') as opfile:
        for itree, smi in enumerate(roots):
            opfile.write(f'\n\n\n--------Tree {itree}: {smi}--------\n')
            
            root_node = TreeNode(smi)
            next_level_nodes = [root_node]
            
            while next_level_nodes:
                for node in next_level_nodes:
                    if node.SMILES in children_dict:
                        children_smis = children_dict[node.SMILES]
                        for child_smi in children_smis:
                            if real_son_node(child_smi, children_smis, parent_dict):
                                node.add_child(TreeNode(child_smi))
                
                tmp_next = []
                for node in next_level_nodes:
                    tmp_next.extend(node.children)
                next_level_nodes = tmp_next
            
            tree_smiles = root_node.print_tree(df_table=df_table, opfile=opfile)
            tree_smile_list.append(tree_smiles)
            tree_list.append(root_node)
    
    # Visualize trees
    logger.info("Generating tree visualizations")
    os.chdir(tree_path)
    Path('Images').mkdir(exist_ok=True, parents=True)
    
    for itree, tree_smiles in enumerate(tree_smile_list):
        d = graphviz.Digraph(filename=root_title)
        d.node_attr["shape"] = "plaintext"
        d.node_attr["fixedsize"] = 'true'
        d.node_attr["height"] = '1'
        d.node_attr["width"] = '1'
        d.node_attr["label"] = ''
        d.attr(fontsize='20')
        d.node_attr.update(fontsize='18')
        d.edge_attr.update(fontsize='18')
        
        level_dummy_nodes = []
        
        for level_idx, level_smiles in enumerate(tree_smiles):
            if level_idx > max_level:
                break
            
            with d.subgraph() as s:
                s.attr(rank='same')
                count = -1
                level_node_dict = {}
                
                for node_idx, cpds_node in enumerate(level_smiles):
                    for cpd in cpds_node:
                        count += 1
                        img_path = f"Images/L{level_idx}_{count}.png"
                        
                        result = show_cpd(
                            [img_path, cpd],
                            df_table=df_table,
                            df_act=df_act,
                            actCols=act_cols,
                            highlightDictList=highlight_dict
                        )
                        
                        if result is None:
                            continue
                        
                        cpd_smi = result[1]
                        
                        if cpd_smi not in level_node_dict:
                            node_label = f"L{level_idx}_{count}"
                            s.node(node_label, image=result[0])
                            
                            if level_idx > 0:
                                d.edge(
                                    f'L{level_idx-1}_{node_idx}',
                                    node_label,
                                    penwidth='0.2',
                                    arrowsize='0.2'
                                )
                            
                            level_node_dict[cpd_smi] = node_label
                        else:
                            node_label = f"L{level_idx}_{count}"
                            s.node(node_label)
                            level_dummy_nodes.append(node_label)
                            
                            existing_label = level_node_dict[cpd_smi]
                            if level_idx > 0 and f'L{level_idx-1}_{node_idx}' not in level_dummy_nodes:
                                d.edge(
                                    f'L{level_idx-1}_{node_idx}',
                                    existing_label,
                                    penwidth='0.2',
                                    arrowsize='0.2'
                                )
        
        d.render(view=False)
    
    logger.info(f"Tree generation complete. Output: {tree_path}")
