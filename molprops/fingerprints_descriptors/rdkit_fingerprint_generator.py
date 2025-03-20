#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

@contents : This module contains the function to generate the Morgan fingerprint from a given SMILES string.
@project :  molprops
@program :  CEU Mass Mediator
@file :  build_data.py
@author :  Alberto Gil De la Fuente (alberto.gilf@gmail.com)
           Pablo Cañadas Miquel (pcanadasmc@gmail.com)
           
@version :  0.0.1, 19 Mar 2025
@information : A valid license of AlvaDesc is necessary to generate the descriptors and fingerprints of chemical structures. 

@copyright :  GNU General Public License v3.0
              Permissions of this strong copyleft license are conditioned on making available complete source code of licensed works and modifications, 
              which include larger works using a licensed work, under the same license. 
              Copyright and license notices must be preserved. Contributors provide an express grant of patent rights.
"""


from typing import Tuple
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolops import FastFindRings


def get_morgan_fingerprint_rdkit(smiles: str):
    """
    Generates the Morgan fingerprint from a given SMILES string.
    
    Syntax:
        fingerprint = get_morgan_fingerprint_rdkit(smiles)
    Parameters:
        smiles (str): A SMILES (Simplified Molecular Input Line Entry System) string representing the chemical structure.
    Returns:
        str: The corresponding Morgan fingerprint of the chemical structure.
    Exceptions:
        ValueError if smiles does not correspond to a molecule.

    """

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    mol.UpdatePropertyCache()
    FastFindRings(mol)

    generator = AllChem.GetMorganGenerator(radius=2, fpSize=1024) 
    fp = generator.GetFingerprint(mol)

    return fp.ToBitString()