#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents : This module contains the exception for when a CAS number cannot be resolved to a structure.
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

class CASNotFoundError(Exception):
    """Exception raised when a CAS number cannot be resolved to a structure."""
    def __init__(self, cas_number, message=None):
        if message is None:
            message = f"CAS number not found: {cas_number}"
        super().__init__(message)
        self.cas_number = cas_number