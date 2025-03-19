#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents :  This module contains the tests for the RDKit molfile downloaders
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
from molprops.molfile_download.molfile_downloader import download_sdf_hmdb, download_sdf_pubchem

class TestMolfileDownloader(unittest.TestCase):

    def test_download_pubchem_sdf(self):
        sdfPath = '/home/ceu/research/repos/cmm_rt_shared/SDF/3D/pubchem'
        pubchem_id = 1
        file_path = os.path.join(self.test_dir, f"{pubchem_id}.sdf")

        expected_fp_length = 1024
        download_sdf_pubchem(1,sdfPath)
        self.assertTrue(os.path.exists(file_path))
        with open(file_path, 'r') as file:
            content = file.read()
            self.assertTrue(len(content) > 0)

