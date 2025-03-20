#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@contents : 
This module contains utility functions for converting data types, specifically for converting a string representation of a large integer into a list of its individual digits.

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

def list_of_ints_from_str(big_int_str: str):
    """
    Convert a string representation of a large integer into a list of its individual digits.

    Args:
        big_int_str (str): A string representing a large integer.

    Returns:
        list: A list of integers, each representing a digit from the input string.
    """
    ints_list = [int(d) for d in str(big_int_str)]
    return ints_list
