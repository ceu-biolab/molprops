#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents : This 
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
import os
import urllib
import time

from molecular_file.sdf_type import SDFType

def download_sdf_pubchem(pc_id, output_path, sdf_type:SDFType = SDFType.THREE_D):
    """ 
        Get SDF file from the pubchem identifier. It retries the call 3 times if the request is not responded.

        Syntax
        ------
          str = download_sdf_pubchem(pc_id, output_path)

        Parameters
        ----------
            [in] pc_id: PC_ID integer corresponding to the pubchem identifier
            [out] output_path: file path to save the corresponding {pc_id}.sdf file
            [in] sdf_type: SDFType enum to specify the type of SDF file [TWO_D, THREE_D]

        Returns
        -------
            None

        Exceptions
        ----------
          Exception:
            If the pubchem identifier is not present in pubchem database it will reraise the exception

        Example
        -------
          >>> inchi_key = get_inchi_key_from_pubchem(1,'.')
    """
    file_path = f"{output_path}/{pc_id}.sdf"
    if os.path.exists(file_path):
        return 
    
    url_pubchem="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + str(pc_id) + "/SDF"
    if sdf_type == SDFType.TWO_D:
        pass
    else:
        url_pubchem = url_pubchem + "?record_type=3d"

    with urllib.request.urlopen(url_pubchem) as response:
        content = response.read().decode("utf-8")

        with open(file_path, "w") as file:
            file.write(content)
            return

def download_sdf_hmdb(hmdb_id, output_path, sdf_type: SDFType = SDFType.THREE_D, retries=3):
    """ 
        Get SDF file from the HMDB identifier. It retries the call 3 times if the request is not responded.

        Syntax
        ------
          None = download_sdf_hmdb(hmdb_id, output_path, sdf_type)

        Parameters
        ----------
            [in] hmdb_id: str
                HMDB identifier (e.g., HMDB0000123).
            [out] output_path: str
                File path to save the corresponding {hmdb_id}.sdf file.
            [in] sdf_type: SDFType (Enum)
                Specifies the type of SDF file [TWO_D, THREE_D].
            [in] retries: int
                Number of retry attempts in case of request failure (default = 3).

        Returns
        -------
            None

        Exceptions
        ----------
          Exception:
            If the HMDB ID is not found in the HMDB database, an exception is raised.

        Example
        -------
          >>> download_sdf_hmdb("HMDB0000123", "./")
    """  

    file_path = f"{output_path}/{hmdb_id}.sdf"
    if os.path.exists(file_path):
        return  # File already exists, no need to download again
    
    # Construct the URL based on SDF type
    url_hmdb = f"https://hmdb.ca/structures/metabolites/{hmdb_id}/download.sdf"
    if sdf_type == SDFType.THREE_D:
        url_hmdb += "?dim=3d"

    # Retry mechanism for network errors
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url_hmdb) as response:
                content = response.read().decode("utf-8")

                with open(file_path, "w") as file:
                    file.write(content)
                    return  # Successfully saved, exit function

        except Exception as e:
            if attempt < retries - 1:
                time.sleep(2)  # Wait before retrying
            else:
                raise Exception(f"Failed to download SDF for HMDB ID {hmdb_id}: {e}")

def main():
    download_sdf_pubchem(1, '.', SDFType.TWO_D)
    download_sdf_hmdb("HMDB0000123", "./", SDFType.THREE_D)

if __name__ == "__main__":
    main()