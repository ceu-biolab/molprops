#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents :  This module contains the tests for the Alvadesc fingerprint generator
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
from molprops.config import ALVADESC_LOCATION
from alvadesccliwrapper.alvadesc import AlvaDesc
from molprops.fingerprints_descriptors.alvadesc_fingerprint_generator import get_fingerprint, get_descriptors
from molprops.fingerprints_descriptors.fingerprint_type import FingerprintType

class TestRDKitFingerprintGenerator(unittest.TestCase):

    def setUp(self):
        aDescPath = ALVADESC_LOCATION
        self.aDesc = AlvaDesc(aDescPath)
        self.sdfPath = '/home/ceu/research/repos/cmm_rt_shared/SDF/3D/pubchem/'
        self.inputFile = "1.sdf"

    def test_get_fingerprint_ecfp(self):
        expResult = "0000001000000000000000000000000000000000000000000000000000000000000000000000010000100000000010000000010000000000001011000000000000010000000000000001000000000000000000000000000000000000000000000000000010000000000000000000000000000000000001000010000001000100000010000000000000001100010000000000000000100000000000000001000010000000000000000000000000000100000100000000000000000000000000010000010000000000000000001001000000000000100000000000000000000000100000001000000001000000000000000010000000010000000000000000000000000000000000000000000100000000000000000100000000000000000000000000000000000100000000000010000000000010001000001001000000010000000000000000000011100000000000000000010000010000000000001000001000000000000000000000000000000000000000100000000000000000000000000000000000000000000001110000000000000000000001000001000000000000000000000000000000010100000000000000000100000000000000000000000000000010000000000000001010000000001000101000000000000000000000100000000000000000000000000000000000010000000000000001000000000000"
        actualResult = get_fingerprint(self.aDesc, self.sdfPath + self.inputFile, fingerprint_type=FingerprintType.ECFP)
        self.assertEqual(expResult, actualResult, "Test FAIL. Check the method get_fingerprint(aDesc, structureFileName, 'ECFP').")

    def test_get_fingerprint_maccs(self):
        expResult = "0000000000000000000000000000010000000000000000001000000000000000000000000100000000001100100010100001000000010001001000000110010000010001000110000101100111101111100100"
        actualResult = get_fingerprint(self.aDesc, self.sdfPath + self.inputFile, fingerprint_type = FingerprintType.MACCSFP)
        self.assertEqual(expResult, actualResult, "Test FAIL. Check the method get_fingerprint(aDesc, structureFileName, 'MACCSFP').")
    
    def test_get_fingerprint_pfp(self):
        expResult = "0000000000000000000000000000000000000000000000000000000000000000000000000000010000100000000010100000010000000000001010000001000000010010000000000000000000010000000000000000000000000100000000000001000000000000000000000000000000000000000000000010100001000000000010000100010000001100010000000000001000100000000100000000100011000000000000100000000001000100000000000000000000000000000000010000010000000000000000001001000000000000100000000000000000000000100000001010000001000000000000000000000101001000000000000000000000000000000000000000000100010000000000000100000000000000000000010010000000000100000000000010000000000010001000001000000000010000000000000000000111000000000000001000000000010000000000000000001000000000000000000000000000000010010000110100001000000000000000000000000000000000000000110000000000000000000001000001000000000000000000000000000000000100000000000000000000100000000010000100000000000010000000010000000010001000001000100000000000000000000000100000000000000001000000000100000000010000000000000001000000000000"
        actualResult = get_fingerprint(self.aDesc, self.sdfPath + self.inputFile, fingerprint_type = FingerprintType.PFP)
        self.assertEqual(expResult, actualResult, "Test FAIL. Check the method get_fingerprint(aDesc, structureFileName, 'PFP').")

    def test_get_desriptors(self):
        expResult_len = 5666
        actualResult = get_descriptors(self.aDesc, self.sdfPath + self.inputFile)
        self.assertEqual(expResult_len, len(actualResult), "Test FAIL. Check the method get_descriptors(aDesc, structureFileName).")

if __name__ == '__main__':
    unittest.main()
