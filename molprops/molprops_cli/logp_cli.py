#!/usr/bin/env python
# -*- coding: utf-8 -*-

 
"""
@contents : This module provides functions for a command line interface to obtain logP from RDKit
@project :  molprops (CEU Mass Mediator project)
@program :  CEU Mass Mediator
@file :  logp_cli.py
@author :  Alberto Gil De la Fuente (alberto.gilf@gmail.com)
           

@version :  0.0.1, 22 Apr 2025
@information : A valid license of AlvaDesc is necessary to generate the descriptors and fingerprints of chemical structures. 

@copyright :  GNU General Public License v3.0
              Permissions of this strong copyleft license are conditioned on making available complete source code of licensed works and modifications, 
              which include larger works using a licensed work, under the same license. 
              Copyright and license notices must be preserved. Contributors provide an express grant of patent rights.
@end
"""

import os
import argparse
import pandas as pd
from molprops.db_endpoints.cas_ids import get_smiles_from_cas_or_compound_name
from molprops.logP.logP import get_logP_from_smiles
from molprops.exceptions.incorrect_mol import InvalidMolecule
from molprops.exceptions.cas_not_found import CASNotFoundError  # if you defined it

def main():
    """
    Main function to generate LOGPs using RDKit from any CSV file.

    Parameters
    ----------
        [in] input_path: str
            Path containing the input file.
        [in] input_file: str
            Name of the input file.
        [in] delimiter: str
            Column delimiter for the input file. Use tab for \t.
        [in] CAS_id_column_name: str
            Name of the CAS ID column.
        [in] output_path: str
            Output path for the result files.

    Returns
    -------
        None

    Exceptions
    ----------
        None
        
    """

    parser = argparse.ArgumentParser(description='Generate LogPs using RDKit from any CSV file.')
    parser.add_argument('--input_path', required=True, help='Input path containing the input file')
    parser.add_argument('--input_file', required=True, help='Name of the input file')
    parser.add_argument('--delimiter', default=',', help='Column delimiter for the input file. Use tab for \t')
    parser.add_argument('--cas_id_column_name', default='CAS', help='Name of the CAS ID column')
    parser.add_argument('--compound_name_id_column_name', default='CAS', help='Name of the compound_name column')
    parser.add_argument('--output_path', required=True, help='Output path for the result files')
    

    args = parser.parse_args()

    inputPath = args.input_path
    inputFileName = os.path.join(inputPath, args.input_file)
    outputPath = args.output_path
    cas_id_column_name = args.cas_id_column_name
    compound_name_id_column_name = args.compound_name_id_column_name
    delimiter = args.delimiter
    if delimiter == 'tab':
        delimiter = '\t'
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)

    df = pd.read_csv(inputFileName, delimiter=delimiter)
    if cas_id_column_name not in df.columns:
        raise ValueError(f"Column '{cas_id_column_name}' not found in input file")
    logp_results = []
    for index, row in df.iterrows():
        cas = row[cas_id_column_name]
        compound_name = row[compound_name_id_column_name]
        try:
            smiles = get_smiles_from_cas_or_compound_name(cas, compound_name=compound_name)
            logp = get_logP_from_smiles(smiles)
        except (CASNotFoundError, InvalidMolecule, ValueError) as e:
            print(f"[WARN] Could not calculate logP for CAS {cas}: {str(e)}")
            logp = None
        logp_results.append(logp)

    df['logP'] = logp_results

    output_file = os.path.join(outputPath, f"logp_{args.input_file}")
    df.to_csv(output_file, index=False, sep=delimiter)
    print(f"[INFO] Results saved to {output_file}")

if __name__ == "__main__":
    main()