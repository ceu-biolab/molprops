#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents :  This module contains the tests for the ClassyFire functions
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
import urllib
from molprops.classyfire.classyfire_wrapper import get_ancestors_from_classyfire, \
is_a_lipid_from_classyfire, get_lm_id_from_inchi_key, is_in_lipidMaps

class TestClassyFireFunctions(unittest.TestCase):

    def test_get_ancestors_from_classyfire_success(self):
        """Test fetching ancestors successfully"""

        # Replace these with actual InChI and InChIKey values
        inchi = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
        inchi_key = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        
        try:
            ancestors = get_ancestors_from_classyfire(inchi, inchi_key)
            self.assertIsInstance(ancestors, list) 
            self.assertGreater(len(ancestors), 0) 
            expected_ancestors = [
            "Acylsalicylic acids", "Acylsalicylic acids and derivatives", "Benzene and substituted derivatives",
            "Benzenoids", "Benzoic acids", "Benzoic acids and derivatives", "Benzoyl derivatives", "Carbonyl compounds",
            "Carboxylic acid derivatives", "Carboxylic acid esters", "Carboxylic acids", "Carboxylic acids and derivatives",
            "Chemical entities", "Dicarboxylic acids and derivatives", "Hydrocarbon derivatives", "Hydroxybenzoic acid derivatives",
            "Organic acids and derivatives", "Organic compounds", "Organic oxides", "Organic oxygen compounds",
            "Organooxygen compounds", "Phenol esters", "Phenoxy compounds", "Salicylic acid and derivatives"
            ]
            for ancestor in expected_ancestors:
                with self.subTest(ancestor=ancestor):
                    self.assertIn(ancestor, ancestors, f"Ancestor '{ancestor}' is missing from the list.")

        except Exception as e:
            self.fail(f"Test failed with exception: {e}")

    def test_get_ancestors_from_classyfire_not_found(self):
        """Test handling of a compound not found in ClassyFire"""

        inchi = "non_existing_inchi"
        inchi_key = "non_existing_inchi_key"    
    
        with self.assertRaises(urllib.error.HTTPError):
            get_ancestors_from_classyfire(inchi, inchi_key)
        

    def test_is_a_lipid_from_classyfire_true(self):
        """Test if a compound is correctly classified as a lipid"""

        inchi = "InChI=1S/C19H38O4/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-19(22)23-17-18(21)16-20/h18,20-21H,2-17H2,1H3" 
        inchi_key = "QHZLMUACJMDIAE-UHFFFAOYSA-N"
        
        try:
            result = is_a_lipid_from_classyfire(inchi, inchi_key)
            self.assertTrue(result)  # The result should be True if it's a lipid
        except Exception as e:
            self.fail(f"Test failed with exception: {e}")

    def test_is_a_lipid_from_classyfire_false(self):

        inchi = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12" 
        inchi_key = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        
        try:
            result = is_a_lipid_from_classyfire(inchi, inchi_key)
            self.assertFalse(result)  # The result should be True if it's a lipid
        except Exception as e:
            self.fail(f"Test failed with exception: {e}")

    def test_is_a_lipid_from_classyfire_not_found(self):
        """Test handling of a compound not found in ClassyFire"""

        # Use an InChIKey that is unlikely to exist in the database
        inchi = "non_existing_inchi"
        inchi_key = "non_existing_inchi_key"
        
        with self.assertRaises(urllib.error.HTTPError):
            is_a_lipid_from_classyfire(inchi, inchi_key)


    def test_get_lm_id_from_inchi_key_success(self):
        # Use a known INCHI Key that is present in LipidMaps
        inchi_key = "RDHQFKQIGNGIED-UHFFFAOYSA-N"
        
        # Call the function
        lm_id = get_lm_id_from_inchi_key(inchi_key)
        
        # Check that the returned lm_id is a string
        self.assertEqual(lm_id, "LMFA07070060", "LM ID should be LMFA07070060")

    def test_get_lm_id_from_inchi_key_failure(self):
        # Use an INCHI Key that is not in LipidMaps
        inchi_key = "INVALIDINCHKEY123"
        
        # Call the function and assert that it raises ValueError
        with self.assertRaises(ValueError):
            get_lm_id_from_inchi_key(inchi_key)

    def test_is_in_lipidMaps_success(self):
        # Use a known INCHI Key that is present in LipidMaps
        inchi_key = "RDHQFKQIGNGIED-UHFFFAOYSA-N"
        
        # Call the function
        result = is_in_lipidMaps(inchi_key)
        
        # Assert that the result is True (INCHI Key exists in LipidMaps)
        self.assertTrue(result)

    def test_is_in_lipidMaps_failure(self):
        # Use an INCHI Key that is not in LipidMaps
        inchi_key = "INVALIDINCHKEY123"
        
        # Call the function
        result = is_in_lipidMaps(inchi_key)
        
        # Assert that the result is False (INCHI Key does not exist in LipidMaps)
        self.assertFalse(result)

if __name__ == "__main__":
    unittest.main()