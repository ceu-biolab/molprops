#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents :  This module contains the tests for the logP calculation. 
@project :  molprops
@program :  CEU Mass Mediator
@file :  classyfire_wrapper.py
@author :  Alberto Gil De la Fuente (alberto.gilf@gmail.com)
           
@version :  0.0.1, 22 Apr 2025
@information : A valid license of AlvaDesc is necessary to generate the descriptors and fingerprints of chemical structures. 

@copyright :  GNU General Public License v3.0
              Permissions of this strong copyleft license are conditioned on making available complete source code of licensed works and modifications, 
              which include larger works using a licensed work, under the same license. 
              Copyright and license notices must be preserved. Contributors provide an express grant of patent rights.
"""

import unittest
from molprops.logP.logP import get_logP_from_smiles
from molprops.exceptions.incorrect_mol import InvalidMolecule

class TestLogPCalculation(unittest.TestCase):

    def test_valid_smiles(self):
        # Oxaloacetic acid: OC(=O)CC(=O)C(O)=O
        smiles = "OC(=O)CC(=O)C(O)=O"
        logp = get_logP_from_smiles(smiles)
        # Check it's a float and within reasonable range
        self.assertIsInstance(logp, float)
        self.assertTrue(-1.0 < logp < -0.5, f"LogP value {logp} is out of expected range")

    def test_invalid_smiles(self):
        # Invalid: gibberish SMILES
        invalid_smiles = "XYZ"
        with self.assertRaises(InvalidMolecule):
            get_logP_from_smiles(invalid_smiles)

if __name__ == "__main__":
    unittest.main()