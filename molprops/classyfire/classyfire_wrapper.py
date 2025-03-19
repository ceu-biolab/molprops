#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

@contents :  This module contains functions to generate fingerprints and descriptors using alvaDesc program
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


import urllib.request,urllib.error, json 
import requests
import time


def get_ancestors_from_classyfire(inchi, inchi_key):
    """
    Get the ancestors from the classyfire classification. 
    It uses the inchi key first, if not classified, inchi.
    Syntax
    -------
            str = get_ancestors_from_classyfire(inchi, inchi_key)
                
    Parameters
    ----------
    [in] inchi: string with the InChI of a compound
    [in] inchi_key: string with the InChIKey of a compound
    
    Returns
    -------
    List of strings representing the ancestors of the compound.
    
    Exceptions
    ----------
    Raises an exception if the request fails or data is unavailable.
    """
    url_classyfire = "http://classyfire.wishartlab.com/entities/" + inchi_key + ".json"
    retries = 0
    while True:
        try:
            with urllib.request.urlopen(url_classyfire) as jsonclassyfire:
                response_code = jsonclassyfire.getcode()

                if response_code == 200:
                    response = requests.get(url_classyfire)
                    response.raise_for_status()  # Raise an HTTPError for bad responses
                    data = response.json()
                    
                    if "ancestors" in data:
                        return data["ancestors"]

        except urllib.error.HTTPError as exception:
            if exception.code == 429:
                if retries < 3:
                    print("Too many requests. Try in 5 seconds")
                    print(exception)
                    time.sleep(5)
            else:
                raise exception


def is_a_lipid_from_classyfire(inchi, inchi_key):
    """ 
        check if the inchi key is a lipid according to classyfire classification. First it uses the gnps2 classyfire endpoint. If it is not there, it goes to the classyfire one. It uses inchi key first, if not classified, inchi. 

        Syntax
        ------
          boolean = is_a_lipid_from_classyfire(inchi, inchi_key)

        Parameters
        ----------

            [in] inchi: string with the inchi of a compound
            [in] inchi_key: string with the inchi key of a compound

        Returns
        -------
          boolean stating wether this inchi key is a lipid according to classyifire classification

        Exceptions
        ----------
          None

        Example
        -------
          >>> is_a_lipid_from_classyfire = is_a_lipid_from_classyfire("RDHQFKQIGNGIED-UHFFFAOYSA-N")
    """

    # http://classyfire.wishartlab.com/entities/BSYNRYMUTXBXSQ-UHFFFAOYSA-N.json
    # https://structure.gnps2.org/classyfire?inchikey=InChIKey=RYYVLZVUVIJVGH-UHFFFAOYSA-N

    url_gnps2_classyfire = "https://structure.gnps2.org/classyfire?inchikey=" + inchi_key
    tried_gnsp2_inchi_key = False
    tried_gnsp2_inchi = False
    retries = 0
    while True:
        try:
            with urllib.request.urlopen(url_gnps2_classyfire) as jsonclassyfire:
                response_code = jsonclassyfire.getcode()

                if response_code == 200:
                    data = json.load(jsonclassyfire)
                    superclass = data["superclass"]["name"]
                    if superclass == "Lipids and lipid-like molecules":
                        return True
                    else:
                        return False
                    #print(data)
                # If the data is not in gnps2, then we hit the classyfire endpoint
                else:
                    if not tried_gnsp2_inchi_key:
                        tried_gnsp2_inchi_key = True
                        url_gnps2_classyfire = "https://structure.gnps2.org/classyfire?inchi=" + inchi
                    elif not tried_gnsp2_inchi:
                        tried_gnsp2_inchi = True
                        url_gnps2_classyfire = "http://classyfire.wishartlab.com/entities/" + inchi_key + ".json"
        except urllib.error.HTTPError as exception:
            if exception.code == 400 or exception.code == 500 or exception.code == 404:
                if not tried_gnsp2_inchi_key:
                    tried_gnsp2_inchi_key = True
                    url_gnps2_classyfire = "https://structure.gnps2.org/classyfire?inchi=" + inchi
                elif not tried_gnsp2_inchi:
                    tried_gnsp2_inchi = True
                    url_gnps2_classyfire = "http://classyfire.wishartlab.com/entities/" + inchi_key + ".json"
                else:
                    raise exception
            elif exception.code == 429:
                if retries < 3:
                    print("Too many requests. Try in 5 seconds")
                    print(exception)
                    time.sleep(5)
            else:
                raise exception


def get_lm_id_from_inchi_key(inchi_key):
    """ 
        Get lm_id from inchi key

        Syntax
        ------
          str = get_lm_id_from_inchi_key(inchi_key)

        Parameters
        ----------
            [in] inchi_key: string with the inchi key of a compound

        Returns
        -------
          str containing the lm id

        Exceptions
        ----------
          Exception:
            ValueError: If the inchi key is not in the lipid maps database
            HttpError: if the url cannot be resolved or the lipidmaps server is down

        Example
        -------
          >>> lm_id = get_lm_id_from_inchi_key("RDHQFKQIGNGIED-UHFFFAOYSA-N")
    """
    url_lipidmaps="https://www.lipidmaps.org/rest/compound/inchi_key/" + inchi_key + "/all"
    while True:
        with urllib.request.urlopen(url_lipidmaps) as jsonLipidMaps:
            data = json.load(jsonLipidMaps)
            if 'lm_id' not in data:
                raise ValueError('LM ID from INCHI KEY ' + inchi_key + ' NOT FOUND')
            else:
                lm_id = data["lm_id"]
                return lm_id

def is_in_lipidMaps(inchi_key):
    """ 
        check if the inchi key is present in lipidmaps

        Syntax
        ------
          boolean = is_in_lipidMaps(inchi_key)

        Parameters
        ----------
            [in] inchi_key: string with the inchi key of a compound

        Returns
        -------
          boolean stating wether the inchi key is present in lipidmaps or not

        Exceptions
        ----------
          Exception if the lipidmaps endpoint is not accesible

        Example
        -------
          >>> is_in_lipidMaps = is_in_lipidMaps("RDHQFKQIGNGIED-UHFFFAOYSA-N")
    """
    try:
        lm_id = get_lm_id_from_inchi_key(inchi_key)
        return True
    except ValueError as ve:
        return False
    except Exception as e:
        raise e