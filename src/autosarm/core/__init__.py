"""
Core modules for AutoSARM.
"""

from autosarm.core.sarm import (
    create_sarm_matrix,
    fragmentize,
    frag_mol_near_ring,
    frags_count,
    attach_left_right,
    has_match,
    sort_mol_sim,
    sort_mol_sim_df,
    key_filter,
    replace_nH,
    match_frag,
    connect_R1,
    get_activity_info,
)

from autosarm.core.tree import (
    create_sar_tree,
    TreeNode,
    find_children,
    find_parents,
    find_root,
    is_parent,
    is_child,
    match_core,
)

__all__ = [
    # SARM functions
    "create_sarm_matrix",
    "fragmentize",
    "frag_mol_near_ring",
    "frags_count",
    "attach_left_right",
    "has_match",
    "sort_mol_sim",
    "sort_mol_sim_df",
    "key_filter",
    "replace_nH",
    "match_frag",
    "connect_R1",
    "get_activity_info",
    # Tree functions
    "create_sar_tree",
    "TreeNode",
    "find_children",
    "find_parents",
    "find_root",
    "is_parent",
    "is_child",
    "match_core",
]
