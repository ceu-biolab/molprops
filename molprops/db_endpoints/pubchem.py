#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents : This module contains the function to get the lipidmaps id from an inchi key
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
import urllib
import json

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
        