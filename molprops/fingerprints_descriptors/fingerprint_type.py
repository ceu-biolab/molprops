#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents :  This module contains the Enum class to specify the type of fingerprint. ECFP, MACCSFP or PFP are from AlvaDesc
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

from enum import Enum

class FingerprintType(Enum):
    '''
    Enum class to specify the type of fingerprint. ECFP, MACCSFP or PFP are from AlvaDesc
    MORGAN are RDKIT fingerprints
    MAP4 are from Reymond group (https://github.com/reymond-group/map4)
    MHFP6 are from Reymond group (https://github.com/reymond-group/mhfp)
    '''
    ECFP = "ECFP"
    MACCSFP = "MACCSFP"
    PFP = "PFP"
    MORGAN = "MORGAN"
    MAP4 = "MAP4"
    MHFP6 = "MHFP6"
