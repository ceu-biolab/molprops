import unittest
from rdkit import Chem
from molprops.molecular_file.rdkit_molfile_writer import get_inchi_from_smiles, inchi_to_3d_structure_rdkit, inchi_to_2d_structure_rdkit, write_mol_file
from molprops.molecular_file.sdf_type import SDFType
import os

class TestRdkitMolfileWriter(unittest.TestCase):

    def test_get_inchi_from_smiles(self):
        smiles = "CCO"
        expected_inchi = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
        self.assertEqual(get_inchi_from_smiles(smiles), expected_inchi)

    def test_inchi_to_3d_structure_rdkit(self):
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        mol = inchi_to_3d_structure_rdkit(inchi=inchi)
        self.assertIsInstance(mol, Chem.Mol)
        self.assertEqual(mol.GetNumAtoms(), 9)  # 3 atoms + 6 hydrogens

    def test_inchi_to_2d_structure_rdkit(self):
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        mol = inchi_to_2d_structure_rdkit(inchi=inchi)
        self.assertIsInstance(mol, Chem.Mol)
        self.assertEqual(mol.GetNumAtoms(), 9)  # 3 atoms + 6 hydrogens

    def test_write_mol_file(self):
        output_path = "test_output"
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        file_path = write_mol_file(output_path, inchi=inchi, sdf_type=SDFType.TWO_D)
        self.assertTrue(os.path.exists(file_path))
        os.remove(file_path)
        os.rmdir(output_path)

if __name__ == '__main__':
    unittest.main()
