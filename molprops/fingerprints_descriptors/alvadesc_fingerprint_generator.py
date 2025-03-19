#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents : 
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

from alvadesccliwrapper.alvadesc import AlvaDesc

from molprops.fingerprints_descriptors.fingerprint_type import FingerprintType
from molprops.molecular_file.file_type import get_file_type

def get_fingerprint(aDesc, mol_structure_path=None, smiles=None, fingerprint_type: FingerprintType = FingerprintType.ECFP, fingerprint_size = 1024):
    """ 
        Generate the the specified type Fingerprint from a molecule structure file

        Syntax
        ------
          str = get_fingerprint(aDesc, mol_structure_path, smiles, fingerprint_type, fingerprint_size)

        Parameters
        ----------
            [in] aDesc (AlvaDesc instance: instance of the aDesc client
            [in] mol_structure_path (str): File name containing the Chemical structure represented by smiles, mol, sdf, mol2 or hin
            [in] smiles (str): SMILES representing the molecule
            [in] fingerprint_type (str): 'ECFP' or 'PFP' or 'MACCSFP'
            [in] fingerprint_size (int): it's not used for MACCS and by default is 1024
        Returns
        -------
          str fingerprint

        Exceptions
        ----------
          TypeError:
            If the mol_structure_path is not smiles, mol, sdf, mol2 or hin
            If the fingerprint_type is not ECFP, PFP or MACCSFP

          RuntimeError:
            If aDesc gets an error calculating the fingerprints

        Example
        -------
          >>> pfp_fingerprint = get_fingerprint((ALVADESC_LOCATION),outputPath + "1.sdf", None, 'PFP')
    """
    if mol_structure_path==None and smiles != None:
        aDesc.set_input_SMILES(smiles)
    elif mol_structure_path != None:
        file_type = get_file_type(mol_structure_path)
        aDesc.set_input_file(mol_structure_path, file_type.name)
    else:
        raise ValueError("SDF or SMILES should be specified")
    if not aDesc.calculate_fingerprint(fingerprint_type.value, fingerprint_size):
        raise RuntimeError('AlvaDesc Error ' + aDesc.get_error())
    else:
        fingerprint = aDesc.get_output()[0]
        return fingerprint
    

def get_descriptors(aDesc, mol_structure_path=None, smiles=None):
    """ 
        Generate all the descriptors from a molecule structure file

        Syntax
        ------
          [obj] = get_descriptors(aDesc, mol_structure_path)

        Parameters
        ----------
            [in] aDesc (AlvaDesc instance): instance of the aDesc client
            [in] mol_structure_path (str): File name containing the Chemical structure represented by smiles, mol, sdf, mol2 or hin
            [in] smiles (str): SMILES representing the molecule
        Returns
        -------
          [obj] descriptors

        Exceptions
        ----------
          TypeError:
            If the mol_structure_path is not smiles, mol, sdf, mol2 or hin

          RuntimeError:
            If aDesc gets an error calculating the descriptors

        Example
        -------
          >>> descriptors = get_descriptors(AlvaDesc(ALVADESC_LOCATION),outputPath + "1.sdf")
    """
    if mol_structure_path==None and smiles != None:
        aDesc.set_input_SMILES(smiles)
    elif mol_structure_path != None:
        file_type = get_file_type(mol_structure_path)
        aDesc.set_input_file(mol_structure_path, file_type.name)
    else:
        raise ValueError("SDF or SMILES should be specified")
    if not aDesc.calculate_descriptors('ALL'):
        raise RuntimeError('AlvaDesc Error ' + aDesc.get_error())
    else:
        descriptors = aDesc.get_output()[0]
        return descriptors