"""
Tests for molecular utility functions.
"""

import pytest
from rdkit import Chem


class TestMolUtils:
    """Tests for autosarm.utils.mol_utils module."""
    
    def test_get_mol_valid_smiles(self):
        """Test get_mol with valid SMILES."""
        from autosarm.utils.mol_utils import get_mol
        
        mol = get_mol("CCO")
        assert mol is not None
        assert mol.GetNumAtoms() == 3
    
    def test_get_mol_invalid_smiles(self):
        """Test get_mol with invalid SMILES."""
        from autosarm.utils.mol_utils import get_mol
        
        mol = get_mol("invalid_smiles")
        assert mol is None
    
    def test_get_mol_none_input(self):
        """Test get_mol with None input."""
        from autosarm.utils.mol_utils import get_mol
        
        mol = get_mol(None)
        assert mol is None
    
    def test_get_mol_mol_input(self):
        """Test get_mol with Mol object input."""
        from autosarm.utils.mol_utils import get_mol
        
        input_mol = Chem.MolFromSmiles("CCO")
        result = get_mol(input_mol)
        assert result is not None
        assert Chem.MolToSmiles(result) == "CCO"
    
    def test_canonic_smiles(self):
        """Test SMILES canonicalization."""
        from autosarm.utils.mol_utils import canonic_smiles
        
        # Different representations of ethanol
        assert canonic_smiles("OCC") == "CCO"
        assert canonic_smiles("C(C)O") == "CCO"
        assert canonic_smiles("CCO") == "CCO"
    
    def test_canonic_smiles_invalid(self):
        """Test canonicalization with invalid input."""
        from autosarm.utils.mol_utils import canonic_smiles
        
        assert canonic_smiles("invalid") is None
    
    def test_is_valid_smiles(self):
        """Test SMILES validation."""
        from autosarm.utils.mol_utils import is_valid_smiles
        
        assert is_valid_smiles("CCO") is True
        assert is_valid_smiles("invalid") is False
        assert is_valid_smiles(None) is False
        assert is_valid_smiles("") is False
        assert is_valid_smiles("123") is False
    
    def test_safe_mol_from_smiles(self):
        """Test safe molecule creation."""
        from autosarm.utils.mol_utils import safe_mol_from_smiles
        
        mol = safe_mol_from_smiles("CCO")
        assert mol is not None
        
        # Should handle numeric input gracefully
        assert safe_mol_from_smiles(123) is None
        assert safe_mol_from_smiles(12.5) is None
    
    def test_compute_fingerprint(self):
        """Test fingerprint computation."""
        from autosarm.utils.mol_utils import compute_fingerprint
        
        fp = compute_fingerprint("CCO")
        assert fp is not None
        
        fp_invalid = compute_fingerprint("invalid")
        assert fp_invalid is None
    
    def test_compute_similarity(self):
        """Test similarity computation."""
        from autosarm.utils.mol_utils import compute_similarity
        
        # Identical molecules should have similarity 1.0
        sim = compute_similarity("CCO", "CCO", mode='smi-smi')
        assert sim == 1.0
        
        # Different molecules should have lower similarity
        sim2 = compute_similarity("CCO", "c1ccccc1", mode='smi-smi')
        assert 0 <= sim2 < 1.0
    
    def test_heavy_atom_count(self):
        """Test heavy atom counting."""
        from autosarm.utils.mol_utils import heavy_atom_count
        
        count = heavy_atom_count("CCO")
        assert count == 3  # 2 carbons + 1 oxygen
    
    def test_remove_dummy(self):
        """Test dummy atom removal."""
        from autosarm.utils.mol_utils import remove_dummy
        
        result = remove_dummy("C(*)C")
        assert "*" not in result


class TestKekule:
    """Tests for Kekulization functions."""
    
    def test_kekulize_smi(self):
        """Test SMILES Kekulization."""
        from autosarm.utils.mol_utils import kekulize_smi
        
        aromatic = "c1ccccc1"
        kekule = kekulize_smi(aromatic)
        
        # Kekulized form should not have lowercase letters
        assert 'c' not in kekule


class TestMapper:
    """Tests for parallel processing mapper."""
    
    def test_mapper_single_job(self):
        """Test mapper with single job."""
        from autosarm.utils.mol_utils import mapper
        
        map_func = mapper(1)
        result = list(map_func(lambda x: x * 2, [1, 2, 3]))
        assert result == [2, 4, 6]
