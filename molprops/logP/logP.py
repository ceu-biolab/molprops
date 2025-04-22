#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents : This module contains the logP calculation using RDKit
@project :  molprops
@program :  CEU Mass Mediator
@file :  classyfire_wrapper.py
@author :  Alberto Gil De la Fuente (alberto.gilf@gmail.com)
           
@version :  0.0.1, 21 Apr 2025
@information : A valid license of AlvaDesc is necessary to generate the descriptors and fingerprints of chemical structures. 

@copyright :  GNU General Public License v3.0
              Permissions of this strong copyleft license are conditioned on making available complete source code of licensed works and modifications, 
              which include larger works using a licensed work, under the same license. 
              Copyright and license notices must be preserved. Contributors provide an express grant of patent rights.
""" 

from rdkit import Chem
from rdkit.Chem import Crippen

from molprops.exceptions.incorrect_mol import InvalidMolecule

def get_logP_from_smiles(smiles):
    """
    Calculate the logP value from a SMILES string. 
    The logP value is a measure of the lipophilicity of a compound, which is important in drug design and development.
    The logP value is calculated using the RDKit library.

    """
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Crippen.MolLogP(mol)
    else:
        raise InvalidMolecule(smiles, "Invalid SMILES string")


