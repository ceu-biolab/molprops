#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents :  This module contains the tests for the external db endpoints
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
import os
from rdkit import Chem
from molprops.db_endpoints.cmm import get_hmdb_id
from molprops.db_endpoints.pubchem import get_lm_id_from_inchi_key
class TestMolfileDownloader(unittest.TestCase):

    def test_get_lm_id_from_inchi_key(self):
        inchi_key = "RDHQFKQIGNGIED-UHFFFAOYSA-N"
        
        lm_id = get_lm_id_from_inchi_key(inchi_key)
        self.assertEqual(lm_id, "LMFA07070060")

    def test_get_lm_id_from_inchi_key_not_found(self):
        inchi_key = "INVALID_INCHI_KEY"
        with self.assertRaises(ValueError):
            get_lm_id_from_inchi_key(inchi_key)

    def test_get_hmdb_from_cmm(self):
        inchi = "InChI=1S/C9H17NO4/c1-7(11)14-8(5-9(12)13)6-10(2,3)4/h8H,5-6H2,1-4H3"
        expected_hmdb_id = "HMDB0240773"
        hmdb_id = get_hmdb_id(inchi)
        self.assertEqual(hmdb_id, expected_hmdb_id)

