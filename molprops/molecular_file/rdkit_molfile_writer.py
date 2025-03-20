#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents : This module provides functions to represent chemical structures as 2D or 3D molecules (depending on the molecule passed) 
            and save them as SDF files.
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

from rdkit import Chem
from rdkit.Chem import AllChem
import os
from molprops.molecular_file.sdf_type import SDFType
def get_inchi_from_smiles(smiles: str):
    """
    Generates the InChI (International Chemical Identifier) from a given SMILES string.

    Syntax:
        inchi = get_inchi_from_smiles(smiles)

    Parameters:
        smiles (str): A SMILES (Simplified Molecular Input Line Entry System) string representing the chemical structure.

    Returns:
        str: The corresponding InChI representation of the chemical structure.

    Exceptions:
        None

    Example:
        >>> get_inchi_from_smiles("CCO")
        'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'

    Note:
        This function relies on the RDKit library for chemical structure manipulation.
    """
    mol = Chem.MolFromSmiles(smiles)

    inchi = Chem.MolToInchi(mol, options='-SNon')
    return inchi

def inchi_to_3d_structure_rdkit(inchi: str=None, smiles: str=None):
    """
    Converts an InChI or SMILES string into a 3D molecular structure using RDKit.

    This function takes either an InChI or SMILES string, converts it into an RDKit 
    molecule object, generates 3D coordinates, and optimizes the structure.

    Parameters:
    ----------
    inchi : str, optional
        The InChI (IUPAC International Chemical Identifier) string of the molecule.
    smiles : str, optional
        The SMILES (Simplified Molecular Input Line Entry System) representation 
        of the molecule.

    Returns:
    -------
    mol : rdkit.Chem.Mol
        An RDKit molecule object with 3D coordinates.

    Raises:
    ------
    ValueError:
        - If neither InChI nor SMILES is provided.
        - If the given InChI or SMILES cannot be converted into an RDKit molecule.

    Notes:
    ------
    - If an InChI is provided, it is converted into an RDKit molecule, and the 
      corresponding SMILES string is generated.
    - If a SMILES string is provided, it is converted into an RDKit molecule, and 
      the corresponding InChI string is generated.
    - Hydrogen atoms are explicitly added before 3D generation.
    - The function first attempts to generate 3D coordinates using `AllChem.ETKDG()`.
      If this fails, it retries using random coordinates.
    - The 3D structure is optimized using the UFF (Universal Force Field) method.

    Example:
    --------
    ```python
    from rdkit import Chem

    inchi_str = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
    mol = inchi_to_3d_structure_rdkit(inchi=inchi_str)
    Chem.MolToMolFile(mol, "output.sdf")
    ```
    """
    # Convert InChI to RDKit molecule object
    if inchi is None and smiles is None:
        raise ValueError("Either InChI or SMILES must be provided")
    elif inchi is not None:
        mol = Chem.MolFromInchi(inchi)
        smiles = Chem.MolToSmiles(mol)
    elif smiles is not None:
        mol = Chem.MolFromSmiles(smiles)
        inchi = get_inchi_from_smiles(smiles)
    
    if mol is None:
        raise ValueError("Invalid InChI string")
    
    # Add hydrogen atoms to the molecule
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if result == -1:
        result = AllChem.EmbedMolecule(mol, useRandomCoords=True)
    
    # Optimize the 3D structure
    AllChem.UFFOptimizeMolecule(mol)
    
    return mol

def inchi_to_2d_structure_rdkit(inchi: str=None, smiles: str=None):
    """
    Converts an InChI or SMILES string into a 2D molecular structure using RDKit.

    Parameters:
    ----------
    inchi : str, optional
        The InChI (IUPAC International Chemical Identifier) string of the molecule.
    smiles : str, optional
        The SMILES (Simplified Molecular Input Line Entry System) representation 
        of the molecule.

    Returns:
    -------
    mol : rdkit.Chem.Mol
        An RDKit molecule object with 2D coordinates.

    Raises:
    ------
    ValueError:
        - If neither InChI nor SMILES is provided.
        - If the given InChI or SMILES cannot be converted into an RDKit molecule.

    Notes:
    ------
    - Uses RDKit’s `Compute2DCoords()` to ensure the molecule is represented in 2D.
    - Hydrogen atoms are explicitly added before coordinate generation.
    """
    if inchi is None and smiles is None:
        raise ValueError("Either InChI or SMILES must be provided")

    if inchi is not None:
        mol = Chem.MolFromInchi(inchi)
    else:
        mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError("Invalid InChI or SMILES string")

    # Add hydrogen atoms
    mol = Chem.AddHs(mol)

    # Generate 2D coordinates
    Chem.rdDepictor.Compute2DCoords(mol)

    return mol

def inchi_to_3d_structure_rdkit(inchi:str =None, smiles: str=None):
    """
    Converts an InChI or SMILES string into a 3D molecular structure using RDKit.

    This function takes either an InChI or SMILES string, converts it into an RDKit 
    molecule object, generates 3D coordinates, and optimizes the structure.

    Parameters:
    ----------
    inchi : str, optional
        The InChI (IUPAC International Chemical Identifier) string of the molecule.
    smiles : str, optional
        The SMILES (Simplified Molecular Input Line Entry System) representation 
        of the molecule.

    Returns:
    -------
    mol : rdkit.Chem.Mol
        An RDKit molecule object with 3D coordinates.

    Raises:
    ------
    ValueError:
        - If neither InChI nor SMILES is provided.
        - If the given InChI or SMILES cannot be converted into an RDKit molecule.

    Notes:
    ------
    - If an InChI is provided, it is converted into an RDKit molecule, and the 
      corresponding SMILES string is generated.
    - If a SMILES string is provided, it is converted into an RDKit molecule, and 
      the corresponding InChI string is generated.
    - Hydrogen atoms are explicitly added before 3D generation.
    - The function first attempts to generate 3D coordinates using `AllChem.ETKDG()`.
      If this fails, it retries using random coordinates.
    - The 3D structure is optimized using the UFF (Universal Force Field) method.

    Example:
    --------
    ```python
    from rdkit import Chem

    inchi_str = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
    mol = inchi_to_3d_structure_rdkit(inchi=inchi_str)
    Chem.MolToMolFile(mol, "output.sdf")
    ```
    """
    # Convert InChI to RDKit molecule object
    if inchi is None and smiles is None:
        raise ValueError("Either InChI or SMILES must be provided")
    elif inchi is not None:
        mol = Chem.MolFromInchi(inchi)
        smiles = Chem.MolToSmiles(mol)
    elif smiles is not None:
        mol = Chem.MolFromSmiles(smiles)
        inchi = get_inchi_from_smiles(smiles)
    
    if mol is None:
        raise ValueError("Invalid InChI string")
    
    # Add hydrogen atoms to the molecule
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if result == -1:
        result = AllChem.EmbedMolecule(mol, useRandomCoords=True)
    
    # Optimize the 3D structure
    AllChem.UFFOptimizeMolecule(mol)
    
    return mol


def write_mol_file(output_path: str, inchi: str=None, inchi_key: str=None, smiles: str=None, sdf_type: SDFType=SDFType.THREE_D):
    """
    Converts an InChI or SMILES string to an RDKit molecule and saves it as an SDF file. 
    If the file is already present in the outputPath, it does not overwrite it.
    If the InChI Key is not provided, it is extracted from the InChI string.


    Parameters:
    ----------
    output_path : str
        The directory where the file should be saved.
    inchi : str, optional
        The InChI (IUPAC International Chemical Identifier) string of the molecule.
    inchi_key : str, optional
        The InChI Key of the molecule.
    smiles : str, optional
        The SMILES (Simplified Molecular Input Line Entry System) representation 
        of the molecule.
    sdf_type : SDFType, optional
        The type of SDF file to create (2D or 3D).
    

    Returns:
    -------
    str
        Full path of the saved SDF file.

    Raises:
    ------
    ValueError:
        - If neither InChI nor SMILES is provided.
    FileNotFoundError:
        - If the specified output directory does not exist.

    Notes:
    ------
    - Calls `get_ND_molecule` to generate a a molecule.
    - Extracts the InChI Key and uses it as the filename.
    - If the directory does not exist, an error is raised.

    Example:
    --------
    ```python
    write_mol_file("output_folder", inchi="InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
    ```
    """
    file_name = f"{inchi_key}.sdf"
    full_path = os.path.join(output_path, file_name)
    if os.path.exists(full_path):
        return full_path 
    
    if sdf_type == SDFType.THREE_D:
        mol = inchi_to_3d_structure_rdkit(inchi=inchi, smiles=smiles)
    else:
        mol = inchi_to_2d_structure_rdkit(inchi=inchi, smiles=smiles)

    # Get the InChI Key
    inchi_key = Chem.InchiToInchiKey(Chem.MolToInchi(mol))

    # Construct the filename based on InChI Key
    file_name = f"{inchi_key}.sdf"
    full_path = os.path.join(output_path, file_name)

    if os.path.exists(full_path):
        return full_path  # Return the existing file path

    # Ensure the output directory exists
    if not os.path.exists(output_path):
        raise FileNotFoundError(f"Output directory '{output_path}' does not exist.")

    # Write molecule to SDF file
    Chem.MolToMolFile(mol, full_path)

    return full_path



