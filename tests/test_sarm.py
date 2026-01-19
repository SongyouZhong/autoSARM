"""
Tests for SARM core functionality.
"""

import pytest
import pandas as pd


class TestFragmentation:
    """Tests for molecular fragmentation."""
    
    def test_frag_mol_near_ring_simple(self):
        """Test fragmentation of a simple molecule."""
        from autosarm.core.sarm import frag_mol_near_ring
        
        # Toluene: methyl attached to benzene
        result = frag_mol_near_ring("Cc1ccccc1")
        
        assert result is not None
        assert len(result) > 0
        
        # Should produce fragment pairs
        for pair in result:
            assert len(pair) == 3  # [core, fragment, original]
    
    def test_frag_mol_near_ring_invalid(self):
        """Test fragmentation with invalid input."""
        from autosarm.core.sarm import frag_mol_near_ring
        
        result = frag_mol_near_ring("invalid_smiles")
        assert result is None
    
    def test_fragmentize(self):
        """Test two-round fragmentation."""
        from autosarm.core.sarm import fragmentize
        
        smiles_list = ["Cc1ccccc1", "CCc1ccccc1", "Cc1ccc(C)cc1"]
        
        round1, round2 = fragmentize(smiles_list, n_jobs=1)
        
        assert isinstance(round1, pd.DataFrame)
        assert isinstance(round2, pd.DataFrame)
        assert 'Key' in round1.columns
        assert 'Value' in round1.columns
        assert 'OrgSmi' in round1.columns


class TestMatching:
    """Tests for substructure matching."""
    
    def test_has_match_true(self):
        """Test matching with valid substructure."""
        from autosarm.core.sarm import has_match
        
        assert has_match("Cc1ccccc1", "c1ccccc1") is True
    
    def test_has_match_false(self):
        """Test matching with no substructure."""
        from autosarm.core.sarm import has_match
        
        assert has_match("CCCC", "c1ccccc1") is False
    
    def test_match_frag(self):
        """Test fragment matching."""
        from autosarm.core.sarm import match_frag
        
        assert match_frag("Cc1ccccc1", "c1ccccc1") == 1
        assert match_frag("CCCC", "c1ccccc1") == 0


class TestKeyFilter:
    """Tests for fragment key filtering."""
    
    def test_key_filter_ring_required(self):
        """Test filter requiring ring."""
        from autosarm.core.sarm import key_filter
        
        # Benzene passes ring filter
        assert key_filter("c1ccccc1", ring=True, num_atoms=3) is True
        
        # Hexane fails ring filter
        assert key_filter("CCCCCC", ring=True, num_atoms=3) is False
    
    def test_key_filter_num_atoms(self):
        """Test filter by atom count."""
        from autosarm.core.sarm import key_filter
        
        # Small molecule fails
        assert key_filter("CC", ring=False, num_atoms=5) is False
        
        # Larger molecule passes
        assert key_filter("CCCCCC", ring=False, num_atoms=5) is True
    
    def test_key_filter_aromatic(self):
        """Test filter requiring aromaticity."""
        from autosarm.core.sarm import key_filter
        
        # Aromatic ring passes
        assert key_filter("c1ccccc1", ring=True, num_atoms=3, aromatic=True) is True
        
        # Non-aromatic ring fails
        assert key_filter("C1CCCCC1", ring=True, num_atoms=3, aromatic=True) is False


class TestConnectR1:
    """Tests for R-group connection."""
    
    def test_connect_R1_simple(self):
        """Test connecting R-group to core."""
        from autosarm.core.sarm import connect_R1
        
        # Core with two dummy atoms
        core = "*c1ccc(*)cc1"
        R = "*C"  # Methyl with dummy
        
        result = connect_R1(R=R, core=core, left_or_right=0)
        
        # Should produce a valid SMILES
        assert result is not None
        from autosarm.utils.mol_utils import get_mol
        assert get_mol(result) is not None


class TestSortMolSim:
    """Tests for similarity-based sorting."""
    
    def test_sort_mol_sim(self):
        """Test similarity sorting."""
        from autosarm.core.sarm import sort_mol_sim
        
        smiles_list = ["c1ccccc1", "Cc1ccccc1", "CCCCCC"]
        
        result = sort_mol_sim(smiles_list)
        
        assert isinstance(result, pd.DataFrame)
        assert len(result) == len(smiles_list)


class TestReplaceNH:
    """Tests for nitrogen state normalization."""
    
    def test_replace_nH(self):
        """Test [nH] replacement."""
        from autosarm.core.sarm import replace_nH
        
        result = replace_nH("[nH]1cccc1")
        assert "[nH]" not in result
        assert "n" in result
