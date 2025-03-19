#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents :  This module contains the Enum class to specify the file Type. A file can be SMILES, MDL, SYBYL, or HYPERCHEM
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

class FileType(Enum):
    SMILES = "SMILES"
    MDL = "MDL"
    SYBYL = "SYBYL"
    HYPERCHEM = "HYPERCHEM"

def get_file_type(mol_structure_path):
    """ 
    Get the file type from a chemical structure file extension

    Parameters
    ----------
    mol_structure_path: str
        Chemical structure file with extension smiles, mol, sdf, mol2, or hin

    Returns
    -------
    FileType: Enum
        The file type based on the file extension.

    Raises
    ------
    TypeError:
        If the file format is not recognized.
    """
    mol_structure_path = mol_structure_path.lower()

    if mol_structure_path.endswith(".smiles"):
        return FileType.SMILES
    elif mol_structure_path.endswith(".mol") or mol_structure_path.endswith(".sdf"):
        return FileType.MDL
    elif mol_structure_path.endswith(".mol2"):
        return FileType.SYBYL
    elif mol_structure_path.endswith(".hin"):
        return FileType.HYPERCHEM
    else:
        raise TypeError("File formats recognized are smiles, mol, sdf, mol2, or hin")