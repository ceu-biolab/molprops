#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents :  This module contains the tests for the CAS IDs functions
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
from unittest.mock import patch, Mock
from molprops.db_endpoints.cas_ids import get_smiles_from_cas_or_compound_name
from molprops.db_endpoints.cas_ids import get_smiles_from_cas_cactus
from molprops.db_endpoints.cas_ids import get_smiles_from_cas_cts
from molprops.db_endpoints.cas_ids import get_smiles_from_cas_pubchem
from molprops.exceptions.cas_not_found import CASNotFoundError

class TestCasIdsFunctions(unittest.TestCase):
    @patch("molprops.db_endpoints.cas_ids.requests.get")
    def test_get_smiles_from_cactus_success(self, mock_get):
        """Test fetching SMILES successfully"""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = "OC(=O)CC(=O)C(O)=O"
        mock_get.return_value = mock_response
        cas_number = "328-42-7"
        expected_smiles = "OC(=O)CC(=O)C(O)=O"
        result = get_smiles_from_cas_cactus(cas_number)
        self.assertEqual(result, expected_smiles)

    
    @patch("molprops.db_endpoints.cas_ids.requests.get")
    def test_invalid_cas_from_cactus(self, mock_get):
        # Setup mock response for invalid CAS
        mock_response = Mock()
        mock_response.status_code = 404
        mock_get.return_value = mock_response

        cas_number = "328-42-71"
        with self.assertRaises(CASNotFoundError):
            get_smiles_from_cas_cactus(cas_number)

    @patch('molprops.db_endpoints.cas_ids.requests.post')
    def test_get_smiles_from_cts_success(self, mock_post):
        mock_post.return_value.status_code = 200
        mock_post.return_value.json.return_value = [{"result": "OC(=O)CC(=O)C(O)=O"}]

        smiles = get_smiles_from_cas_cts("328-42-7")
        self.assertEqual(smiles, "OC(=O)CC(=O)C(O)=O")

    @patch('molprops.db_endpoints.cas_ids.requests.post')
    def test_get_smiles_from_cts_not_found(self, mock_post):
        mock_post.return_value.status_code = 404

        with self.assertRaises(CASNotFoundError):
            get_smiles_from_cas_cts("9999-99-9")

    @patch('molprops.db_endpoints.cas_ids.requests.get')
    def test_get_smiles_from_pubchem_success(self, mock_get):
        mock_response = {
            "PropertyTable": {
                "Properties": [
                    {"CID": 1234, "CanonicalSMILES": "OC(=O)CC(=O)C(O)=O"}
                ]
            }
        }
        mock_get.return_value.status_code = 200
        mock_get.return_value.json.return_value = mock_response

        smiles = get_smiles_from_cas_pubchem("328-42-7")
        self.assertEqual(smiles, "OC(=O)CC(=O)C(O)=O")

    @patch('molprops.db_endpoints.cas_ids.requests.get')
    def test_get_smiles_from_pubchem_not_found(self, mock_get):
        mock_get.return_value.status_code = 404

        with self.assertRaises(CASNotFoundError):
            get_smiles_from_cas_pubchem("328-42-71")
    
    def test_invalid_cas(self):
        cas_number = "328-42-71"  # Invalid CAS
        with self.assertRaises(CASNotFoundError):
            get_smiles_from_cas_or_compound_name(cas_number)

    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_cas_cactus")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_cas_pubchem")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_cas_cts")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_name_pubchem")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_name_cts")
    def test_get_smiles_from_cas_success_cactus(self, mock_name_cts, mock_name_pubchem, mock_cts, mock_pubchem, mock_cactus):
        """Test successful retrieval of SMILES from Cactus"""
        mock_cactus.return_value = "OC(=O)CC(=O)C(O)=O"
        cas_number = "328-42-7"
        result = get_smiles_from_cas_or_compound_name(cas_number)
        self.assertEqual(result, "OC(=O)CC(=O)C(O)=O")
        mock_cactus.assert_called_once_with(cas_number)
        mock_pubchem.assert_not_called()
        mock_cts.assert_not_called()
        mock_name_pubchem.assert_not_called()
        mock_name_cts.assert_not_called()

    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_cas_cactus")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_cas_pubchem")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_cas_cts")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_name_pubchem")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_name_cts")
    def test_get_smiles_from_cas_success_pubchem(self, mock_name_cts, mock_name_pubchem, mock_cts, mock_pubchem, mock_cactus):
        """Test successful retrieval of SMILES from PubChem"""
        mock_cactus.side_effect = CASNotFoundError("328-42-7", "Not found")
        mock_pubchem.return_value = "OC(=O)CC(=O)C(O)=O"
        cas_number = "328-42-7"
        result = get_smiles_from_cas_or_compound_name(cas_number)
        self.assertEqual(result, "OC(=O)CC(=O)C(O)=O")
        mock_cactus.assert_called_once_with(cas_number)
        mock_pubchem.assert_called_once_with(cas_number, None)
        mock_cts.assert_not_called()
        mock_name_pubchem.assert_not_called()
        mock_name_cts.assert_not_called()

    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_cas_cactus")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_cas_pubchem")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_cas_cts")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_name_pubchem")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_name_cts")
    def test_get_smiles_from_cas_success_cts(self, mock_name_cts, mock_name_pubchem, mock_cts, mock_pubchem, mock_cactus):
        """Test successful retrieval of SMILES from CTS"""
        mock_cactus.side_effect = CASNotFoundError("328-42-7", "Not found")
        mock_pubchem.side_effect = CASNotFoundError("328-42-7", "Not found")
        mock_cts.return_value = "OC(=O)CC(=O)C(O)=O"
        cas_number = "328-42-7"
        result = get_smiles_from_cas_or_compound_name(cas_number)
        self.assertEqual(result, "OC(=O)CC(=O)C(O)=O")
        mock_cactus.assert_called_once_with(cas_number)
        mock_pubchem.assert_called_once_with(cas_number, None)
        mock_cts.assert_called_once_with(cas_number)
        mock_name_pubchem.assert_not_called()
        mock_name_cts.assert_not_called()

    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_cas_cactus")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_cas_pubchem")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_cas_cts")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_name_pubchem")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_name_cts")
    def test_get_smiles_from_cas_success_name_pubchem(self, mock_name_cts, mock_name_pubchem, mock_cts, mock_pubchem, mock_cactus):
        """Test successful retrieval of SMILES from PubChem using compound name"""
        mock_cactus.side_effect = CASNotFoundError("328-42-7", "Not found")
        mock_pubchem.side_effect = CASNotFoundError("328-42-7", "Not found")
        mock_cts.side_effect = CASNotFoundError("328-42-7", "Not found")
        mock_name_pubchem.return_value = "OC(=O)CC(=O)C(O)=O"
        cas_number = "328-42-7"
        compound_name = "Oxaloacetic acid"
        result = get_smiles_from_cas_or_compound_name(cas_number, compound_name)
        self.assertEqual(result, "OC(=O)CC(=O)C(O)=O")
        mock_cactus.assert_called_once_with(cas_number)
        mock_pubchem.assert_called_once_with(cas_number, compound_name)
        mock_cts.assert_called_once_with(cas_number)
        mock_name_pubchem.assert_called_once_with(compound_name)
        mock_name_cts.assert_not_called()

    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_cas_cactus")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_cas_pubchem")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_cas_cts")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_name_pubchem")
    @patch("molprops.db_endpoints.cas_ids.get_smiles_from_name_cts")
    def test_get_smiles_from_cas_not_found(self, mock_name_cts, mock_name_pubchem, mock_cts, mock_pubchem, mock_cactus):
        """Test failure to retrieve SMILES from all sources"""
        mock_cactus.side_effect = CASNotFoundError("328-42-7", "Not found")
        mock_pubchem.side_effect = CASNotFoundError("328-42-7", "Not found")
        mock_cts.side_effect = CASNotFoundError("328-42-7", "Not found")
        mock_name_pubchem.side_effect = ValueError("Not found")
        mock_name_cts.side_effect = ValueError("Not found")
        cas_number = "328-42-7"
        compound_name = "Oxaloacetic acid"
        with self.assertRaises(CASNotFoundError):
            get_smiles_from_cas_or_compound_name(cas_number, compound_name)
        mock_cactus.assert_called_once_with(cas_number)
        mock_pubchem.assert_called_once_with(cas_number, compound_name)
        mock_cts.assert_called_once_with(cas_number)
        mock_name_pubchem.assert_called_once_with(compound_name)
        mock_name_cts.assert_called_once_with(compound_name)

if __name__ == "__main__":
    unittest.main()