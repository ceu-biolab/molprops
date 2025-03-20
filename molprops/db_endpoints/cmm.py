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
def get_hmdb_id(inchi: str):
    """ 
        Retrieve the HMDB ID from a MySQL database using an InChI string.

        Syntax
        ------
          str = get_hmdb_id(inchi)

        Parameters
        ----------
            [in] inchi: str
                InChI string of a compound.

        Returns
        -------
            hmdb_id: str
                Corresponding HMDB ID.

        Exceptions
        ----------
          Exception:
            If the InChI is not found in the database.
            If there is a MySQL connection error.

        Example
        -------
          >>> inchi = "InChI=1S/C9H8O4/c10-6-3-1-2-5(4-6)9(13)12-8-7(11)9/h1-4,7-8,11H"
          >>> hmdb_id = get_hmdb_id(inchi)
          >>> print(hmdb_id)
    """  

    import mysql.connector
    
    try:
        
        # Load credentials
        creds = load_db_credentials()

        # Establish MySQL connection using loaded credentials
        conn = mysql.connector.connect(
            host=creds["host"],
            user=creds["user"],
            password=creds["password"],
            database=creds["database"]
        )
        cursor = conn.cursor()
        
        # SQL Query
        query = """
        SELECT hmdb_id 
        FROM compounds_hmdb ch 
        INNER JOIN compound_identifiers ci 
        ON ci.compound_id = ch.compound_id 
        WHERE inchi LIKE %s;
        """
        cursor.execute(query, (inchi,))  # Execute query with parameter

        # Fetch result
        result = cursor.fetchone()

        # Close connection
        cursor.close()
        conn.close()

        # Check if result is found
        if result:
            return result[0]  # Return HMDB ID
        else:
            raise Exception("HMDB ID not found for the given InChI.")

    except mysql.connector.Error as e:
        raise Exception(f"MySQL Error: {str(e)}")


def load_db_credentials(file_path="molprops/db.credentials"):
    """ 
        Load database credentials from a text file.

        Syntax
        ------
          dict = load_db_credentials(file_path)

        Parameters
        ----------
            [in] file_path: str
                Path to the credentials file.

        Returns
        -------
            credentials: dict
                Dictionary containing database connection details.

        Exceptions
        ----------
          Exception:
            If the file is missing or improperly formatted.

        Example
        -------
          >>> creds = load_db_credentials("db_credentials.txt")
    """  
    credentials = {}
    with open(file_path, "r") as file:
        for line in file:
            key, value = line.strip().split("=")
            credentials[key] = value
    return credentials