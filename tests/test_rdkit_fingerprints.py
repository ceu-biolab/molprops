#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents :  This module contains the tests for the RDKit fingerprint generator
@project :  molprops
@program :  CEU Mass Mediator
@file :  classyfire_wrapper.py
@author :  Alberto Gil De la Fuente (alberto.gilf@gmail.com)
           Pablo Cañadas Miquel (pcanadasmc@gmail.com)
           
@version :  0.0.1, 19 Mar 2025
@information : A valid license of AlvaDesc is necessary to generate the descriptors and fingerprints of chemical structures. 

@copyright :  GNU General Public License v3.0
              Permissions of this strong copyleft license are conditioned on making available complete source code of licensed works and modifications, 
              which include larger works using a licensed work, under the same license. 
              Copyright and license notices must be preserved. Contributors provide an express grant of patent rights.
"""


import unittest
from rdkit import Chem
from molprops.fingerprints_descriptors.rdkit_fingerprint_generator import get_morgan_fingerprint_rdkit

class TestRDKitFingerprintGenerator(unittest.TestCase):

    def test_valid_smiles(self):
        smiles = "CCO"
        expected_fp_length = 1024
        fingerprint = get_morgan_fingerprint_rdkit(smiles)
        self.assertEqual(len(fingerprint), expected_fp_length)
        self.assertTrue(isinstance(fingerprint, str))

    def test_invalid_smiles(self):
        smiles = "invalid_smiles"
        with self.assertRaises(ValueError):
            get_morgan_fingerprint_rdkit(smiles)

    def test_empty_smiles(self):
        smiles = ""
        
        expected_fp_length = 1024
        fingerprint = get_morgan_fingerprint_rdkit(smiles)
        expected_fingerprint = "0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
        self.assertEqual(fingerprint, expected_fingerprint)


if __name__ == '__main__':
    unittest.main()
