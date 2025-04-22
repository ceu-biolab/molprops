#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents : 
This module contains a function to retrieve the SMILES representation of a chemical compound using its CAS number.

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
import requests
from molprops.exceptions.cas_not_found import CASNotFoundError


def get_smiles_from_cas_or_compound_name(cas_number, compound_name=None):
    """
    Retrieve the SMILES (Simplified Molecular Input Line Entry System) representation of a compound 
    using its CAS (Chemical Abstracts Service) number. The function attempts to retrieve the SMILES 
    from multiple sources in a prioritized order.

    Parameters:
        cas_number (str): The CAS number of the compound to retrieve the SMILES for.
        compound_name (str, optional): The name of the compound, used as a fallback if the CAS number 
            lookup fails. Defaults to None.

    Returns:
        str: The SMILES representation of the compound.

    Raises:
        CASNotFoundError: If the SMILES cannot be found using the CAS number or compound name 
            across all available sources.
    """
    try:
        # First, try to get the SMILES from Cactus
        return get_smiles_from_cas_cactus(cas_number)
    except CASNotFoundError:
        # If Cactus fails, try to get the SMILES from PubChem
        try: 
            return get_smiles_from_cas_pubchem(cas_number, compound_name)
        except CASNotFoundError:
            # If PubChem fails, try to get the SMILES from CTS
            try:
                return get_smiles_from_cas_cts(cas_number)
            except CASNotFoundError:
                # If CTS fails, try to get the SMILES from PubChem using the compound name
                try:
                    return get_smiles_from_name_pubchem(compound_name)
                except Exception:
                    # If Pubchem compound name fails, try to get the SMILES from CTS using the compound name
                    try:
                        return get_smiles_from_name_cts(compound_name)
                    except (ValueError, ConnectionError):
                        raise CASNotFoundError(cas_number, f"SMILES not found for CAS {cas_number} or name {compound_name}")
                
def get_smiles_from_cas_cactus(cas_number):
    """
    Retrieve the SMILES (Simplified Molecular Input Line Entry System) representation of a chemical compound
    using its CAS (Chemical Abstracts Service) number from the CACTUS Chemical Identifier Resolver.
    Args:
        cas_number (str): The CAS number of the chemical compound to query.
    Returns:
        str: The SMILES representation of the chemical compound.
    Raises:
        CASNotFoundError: If the CAS number is not found or is invalid.
    """

    url = f"https://cactus.nci.nih.gov/chemical/structure/{cas_number}/smiles"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    else:
        raise CASNotFoundError(cas_number, f"CAS number not found or invalid: {cas_number}")


def get_smiles_from_cas_pubchem(cas_number, compound_name=None):
    """
    Retrieve the Canonical SMILES string for a compound using its CAS number from the PubChem database.
    Args:
        cas_number (str): The CAS number of the compound to look up.
        compound_name (str, optional): The name of the compound (not used in the current implementation).
    Returns:
        str: The Canonical SMILES string of the compound.
    Raises:
        CASNotFoundError: If the CAS number is not found in the PubChem database or if the response
                        does not contain the expected data structure.
    Notes:
        - This function uses the PubChem PUG REST API to fetch the Canonical SMILES.
        - Ensure that the `requests` library is installed and accessible in your environment.
    """
    
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_number}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    
    if response.status_code == 200:
        try:
            data = response.json()
            return data['PropertyTable']['Properties'][0]['CanonicalSMILES']
        except (KeyError, IndexError):
            raise CASNotFoundError(cas_number, f"CAS {cas_number} not found in PubChem")
    raise CASNotFoundError(cas_number, f"CAS {cas_number} not found in PubChem")


def get_smiles_from_name_pubchem(name):
    """
    Query PubChem using a compound name and retrieve the canonical SMILES string.

    Parameters
    ----------
    name : str
        Common or IUPAC name of the compound (e.g., "Oxaloacetic acid").

    Returns
    -------
    str
        Canonical SMILES string from PubChem.

    Raises
    ------
    ValueError
        If no SMILES is found in the PubChem response.
    ConnectionError
        If the API request fails with a non-200 status code.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)

    if response.status_code == 200:
        try:
            data = response.json()
            return data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
        except (KeyError, IndexError):
            raise ValueError(f"SMILES not found for name: {name}")
    elif response.status_code == 404:
        raise ValueError(f"Name '{name}' not found in PubChem")
    else:
        raise ConnectionError(f"PubChem request failed with status code {response.status_code}")

def get_smiles_from_cas_cts(cas_number):
    """
    Try the Chemical Translation Service (CTS) to convert a CAS number to a SMILES string.
    This function sends a POST request to the CTS API to retrieve the corresponding 
    SMILES (Simplified Molecular Input Line Entry System) representation for a given 
    CAS (Chemical Abstracts Service) number.
    Args:
        cas_number (str): The CAS number to be converted to a SMILES string.
    Returns:
        str: The SMILES string corresponding to the given CAS number.
    Raises:
        CASNotFoundError: If the CAS number is not found in the CTS database, or if 
                          the response does not contain a valid SMILES string.
        CASNotFoundError: If the CTS request fails with a status code other than 200 or 404.
    """

    url = "https://cts.fiehnlab.ucdavis.edu/rest/convert/CAS/SMILES"
    response = requests.post(url, json=[cas_number])

    if response.status_code == 200:
        try:
            result = response.json()
            return result[0]['result']  # could be a list of SMILES
        except (KeyError, IndexError):
            raise CASNotFoundError(cas_number, f"SMILES not found in CTS response for CAS {cas_number}")
    elif response.status_code == 404:
        
        raise CASNotFoundError(cas_number, f"CAS {cas_number} not found in CTS")
    else:
        raise CASNotFoundError(cas_number, f"CTS request failed with code {response.status_code}")
    
def get_smiles_from_name_cts(compound_name):
    """
    Retrieve the SMILES (Simplified Molecular Input Line Entry System) representation 
    of a chemical compound given its name using the CTS (Chemical Translation Service) API.
    Args:
        compound_name (str): The name of the chemical compound.
    Returns:
        str: The SMILES representation of the compound.
    Raises:
        ValueError: If the SMILES representation is not found for the given compound name.
        ConnectionError: If the request to the CTS API fails with a non-200 status code.
    """
    url = "https://cts.fiehnlab.ucdavis.edu/rest/convert/Name/SMILES"
    response = requests.post(url, json=[compound_name])
    if response.status_code == 200:
        results = response.json()
        if results and 'result' in results[0]:
            return results[0]['result']
        else:
            raise ValueError(f"SMILES not found for name: {compound_name}")
    else:
        raise ConnectionError(f"CTS request failed with status code {response.status_code}")