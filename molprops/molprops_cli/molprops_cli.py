#!/usr/bin/env python
# -*- coding: utf-8 -*-

 
"""
@contents : This module provides functions for a command line interface to generate fingerprints and descriptors for chemical structures using the AlvaDesc software.
@project :  molprops (CEU Mass Mediator Retention Time)
@program :  CEU Mass Mediator
@file :  build_data_descriptors.py
@author :  Alberto Gil De la Fuente (alberto.gilf@gmail.com)
           

@version :  0.0.1, 20 Mar 2025
@information : A valid license of AlvaDesc is necessary to generate the descriptors and fingerprints of chemical structures. 

@copyright :  GNU General Public License v3.0
              Permissions of this strong copyleft license are conditioned on making available complete source code of licensed works and modifications, 
              which include larger works using a licensed work, under the same license. 
              Copyright and license notices must be preserved. Contributors provide an express grant of patent rights.
@end
"""


from alvadesccliwrapper.alvadesc import AlvaDesc
import csv
import os
from rdkit import Chem
import argparse

from molprops.molecular_file.sdf_type import SDFType
from molprops.molecular_file.file_type import get_file_type, FileType
from molprops.fingerprints_descriptors.fingerprint_type import FingerprintType
from molprops.molfile_download.molfile_downloader import download_sdf_pubchem, download_sdf_hmdb
from molprops.db_endpoints.cmm import get_hmdb_id
from molprops.molecular_file.rdkit_molfile_writer import write_mol_file
from molprops.fingerprints_descriptors.alvadesc_fingerprint_generator import generate_vector_fingerprints, get_descriptors, get_fingerprint
from molprops.fingerprints_descriptors.rdkit_fingerprint_generator import get_morgan_fingerprint_rdkit

ALVADESC_LOCATION = '/usr/bin/alvaDescCLI'

def main():
    """
    Main function to generate fingerprints and descriptors using AlvaDesc from any CSV
    file.

    Parameters
    ----------
        [in] input_path: str
            Path containing the input file.
        [in] input_file: str
            Name of the input file.
        [in] delimiter: str
            Column delimiter for the input file. Use tab for \t.
        [in] smiles_column_name: str
            Name of the SMILES column.
        [in] inchi_column_name: str
            Name of the InChI column.
        [in] pubchem_id_column_name: str
            Name of the pubchem ID column.
        [in] hmdb_id_column_name: str
            Name of the HMDB ID column.
        [in] output_path: str
            Output path for the result files.
        [in] input_file_type: str
            Descriptor type: TWO_D or THREE_D.

    Returns
    -------
        None

    Exceptions
    ----------
        None
        
    """

    parser = argparse.ArgumentParser(description='Generate fingerprints and descriptors using AlvaDesc from any CSV file.')
    parser.add_argument('--input_path', required=True, help='Input path containing the input file')
    parser.add_argument('--input_file', required=True, help='Name of the input file')
    parser.add_argument('--delimiter', default=',', help='Column delimiter for the input file. Use tab for \t')
    parser.add_argument('--smiles_column_name', default=None, help='Name of the SMILES column')
    parser.add_argument('--inchi_column_name', default=None, help='Name of the InChI column')
    parser.add_argument('--pubchem_id_column_name', default=None, help='Name of the pubchem ID column')
    parser.add_argument('--hmdb_id_column_name', default=None, help='Name of the HMDB ID column')
    parser.add_argument('--output_path', required=True, help='Output path for the result files')
    parser.add_argument('--input_file_type', default='THREE_D', choices=['TWO_D', 'THREE_D'], help='Descriptor type: TWO_D or THREE_D')
    

    args = parser.parse_args()

    inputPath = args.input_path
    inputFileName = os.path.join(inputPath, args.input_file)
    outputPath = args.output_path
    smiles_column_name = args.smiles_column_name
    inchi_column_name = args.inchi_column_name
    pubchem_id_column_name = args.pubchem_id_column_name
    hmdb_id_column_name = args.hmdb_id_column_name
    input_file_type = args.input_file_type
    input_file_type = getattr(SDFType, input_file_type)
    delimiter = args.delimiter
    if delimiter == 'tab':
        delimiter = '\t'
    sdf_path = "/home/ceu/research/repos/cmm_rt_shared/SDF"
    sdf_path_2D = os.path.join(sdf_path, "2D")
    os.makedirs(sdf_path_2D, exist_ok=True)
    sdf_path_3D = os.path.join(sdf_path, "3D")
    os.makedirs(sdf_path_3D, exist_ok=True)

    if input_file_type == SDFType.TWO_D:
        sdf_path = sdf_path_2D
        outputPath = os.path.join(outputPath, "2D")
    else:
        sdf_path = sdf_path_3D
        outputPath = os.path.join(outputPath, "3D")
    
    
    vector_fingerprints_path = os.path.join(outputPath, "vector_fingerprints")
    os.makedirs(vector_fingerprints_path, exist_ok=True)
    
    #Constants
    NUMBER_FPVALUES=2214
    # VARIABLES OF AlvaDesc Software
    aDesc = AlvaDesc(ALVADESC_LOCATION)

    column_name_is_3D = "is_3D"
    
    # IT WILL TAKE SMILES to create a CSV file containing the vector with fingerprints (ECFP, MACCSFP and PFP) of each compound
    base_fileName, fileExtension = os.path.splitext(args.input_file)
    
    outputFileDescriptorsName = os.path.join(os.path.join(vector_fingerprints_path, base_fileName + "_descriptors" + fileExtension))
    outputFileDescriptorsAndFingerprintsName = os.path.join(os.path.join(vector_fingerprints_path, base_fileName + "_descriptorsAndFingerprints" + fileExtension))
    outputFileFingerprintsVectorizedName = os.path.join(os.path.join(vector_fingerprints_path, base_fileName + "_vectorfingerprintsVectorized" + fileExtension))
    if os.path.isfile(outputFileDescriptorsName):
        os.remove(outputFileDescriptorsName)
    if os.path.isfile(outputFileDescriptorsAndFingerprintsName):
        os.remove(outputFileDescriptorsAndFingerprintsName)
    if os.path.isfile(outputFileFingerprintsVectorizedName):
        os.remove(outputFileFingerprintsVectorizedName)

    with open(inputFileName) as csvfile:
        reader = csv.DictReader(csvfile,delimiter=delimiter,quotechar='"')
        
        # RUN A MOCK SDF TO OBTAIN DESCRIPTORS HEADERS
        inchi="O=C(NCc1ccc(cc1)F)NCCCN1CCc2c1cccc2"
        aDesc.set_input_SMILES(inchi)
        aDesc.calculate_descriptors('ALL')
        listDescriptors = aDesc.get_output_descriptors()

        # Create here the headers from the input file and then add the descriptors
        descriptorFieldNames =reader.fieldnames.copy()
        descriptorFieldNames.append(column_name_is_3D)
        descriptorFieldNames.extend(listDescriptors)
        descriptorsAndFingerPrintsFieldNames = descriptorFieldNames[:]
        
        # Write headers in the output file
        outputFileDescriptors = open(outputFileDescriptorsName, 'w', newline='')
        writerDescriptors = csv.DictWriter(outputFileDescriptors, fieldnames = descriptorFieldNames)
        writerDescriptors.writeheader()

        
        # WRITER FOR FINGERPRINTS AND DESCRIPTORS
        descriptorsAndFingerPrintsFieldNames.append('ECFP')
        descriptorsAndFingerPrintsFieldNames.append('MACCSFP')
        descriptorsAndFingerPrintsFieldNames.append('PFP')
        descriptorsAndFingerPrintsFieldNames.append('MorganFP')

        outputFileDescriptorsAndFingerprints = open(outputFileDescriptorsAndFingerprintsName, 'w', newline='')
        writerDescriptorsAndFingerprints = csv.DictWriter(outputFileDescriptorsAndFingerprints, fieldnames = descriptorsAndFingerPrintsFieldNames)
        writerDescriptorsAndFingerprints.writeheader()
        
        # Create here the headers from the input file and then add the Fingerprints
        FPVectorizedFieldNames = reader.fieldnames.copy()
        FPVectorizedFieldNames.append(column_name_is_3D)
        for i in range(0,NUMBER_FPVALUES):
            header_name = "V" + str(i+1)
            FPVectorizedFieldNames.append(header_name)

        # WRITER FOR FINGERPRINTS VECTORIZED
        outputFileFingerprintsVectorized = open(outputFileFingerprintsVectorizedName, 'w', newline='')
        writerFingerprintsVectorized = csv.DictWriter(outputFileFingerprintsVectorized, fieldnames = FPVectorizedFieldNames)
        writerFingerprintsVectorized.writeheader()

        
        
        descriptors_dict = {}
        maccsfp_dict = {}
        ecfp_dict = {}
        pfp_dict = {}
        morganfp_dict = {}
        vector_fingerprints_dict = {}

        for row in reader:
            pc_id = None
            if pubchem_id_column_name:
                pc_id = row[pubchem_id_column_name]
            
            if pc_id:
                try:
                    pc_id_sdf_path = f"{sdf_path}/pubchem/{pc_id}.sdf"
                    if not os.path.exists(pc_id_sdf_path):
                        download_sdf_pubchem(pc_id,pc_id_sdf_path, sdf_type=input_file_type)
                    
                except Exception as e:
                    pc_id_sdf_path = None
            else:
                pc_id_sdf_path = None
            
            hmdb_id = None
            if hmdb_id_column_name:
                hmdb_id = row[hmdb_id_column_name]
                if not hmdb_id:
                    try:
                        hmdb_id = get_hmdb_id(inchi)
                    except Exception as e:
                        hmdb_id = None
            if hmdb_id:
                try:
                    hmdb_id_sdf_path = f"{sdf_path}/hmdb/{hmdb_id}.sdf"
                    if not os.path.exists(hmdb_id_sdf_path):
                        download_sdf_hmdb(hmdb_id,hmdb_id_sdf_path, sdf_type=input_file_type)
                except Exception as e:
                    hmdb_id_sdf_path = None
            else:
                hmdb_id_sdf_path = None

            if inchi_column_name:
                inchi = row[inchi_column_name]
                if inchi != None:
                    mol = Chem.MolFromInchi(inchi)
                    if not mol:
                        continue
                    smiles = Chem.MolToSmiles(mol)
                    
            elif smiles_column_name:
                smiles = row[smiles_column_name]
                if not smiles:
                    continue
                mol = Chem.MolFromSmiles(smiles)
                if not mol:
                    continue
                inchi = Chem.MolToInchi(mol)
            
            inchi_key = Chem.MolToInchiKey(mol)
            
            if input_file_type == SDFType.THREE_D:
                try:
                    sdf_full_path = write_mol_file(sdf_path_3D, inchi=inchi, inchi_key=inchi_key, smiles=smiles, sdf_type=SDFType.THREE_D)
                    is_3D = True
                    row[column_name_is_3D] = True
                except ValueError as e:
                    sdf_full_path = write_mol_file(sdf_path_2D, inchi=inchi, inchi_key=inchi_key, smiles=smiles, sdf_type=SDFType.TWO_D))
                    is_3D = False
                    row[column_name_is_3D] = False
            else: 
                sdf_full_path = write_mol_file(sdf_path_2D, inchi=inchi, inchi_key=inchi_key, smiles=smiles, sdf_type=SDFType.TWO_D))
                is_3D = False
                row[column_name_is_3D] = False
            
            # Do directly the copy of all elements of the row
            if inchi_key in descriptors_dict:
                descriptors = descriptors_dict[inchi_key]
            else:
                descriptors = get_descriptors(aDesc, mol_structure_path=sdf_full_path, smiles=smiles)
                descriptors_dict[inchi_key] = descriptors
            partialDictDescriptorsRow = row.copy()
            
            for i in range(0,len(listDescriptors)):
                descriptor_header = listDescriptors[i]
                partialDictDescriptorsRow[descriptor_header] = descriptors[i]
            writerDescriptors.writerow(partialDictDescriptorsRow)
            partialDictDescriptorsAndFingerprintsRow = partialDictDescriptorsRow.copy()
            
            # Add fingerprints
            if inchi_key in ecfp_dict:
                fingerprint_ecfp = ecfp_dict[inchi_key]
                fingerprint_maccs = maccsfp_dict[inchi_key]
                fingerprint_pfp = pfp_dict[inchi_key]
                fingerprint_morgan = morganfp_dict[inchi_key]
                vector_fingerprints = vector_fingerprints_dict[inchi_key]
            else:
                
                fingerprint_ecfp = get_fingerprint(aDesc,mol_structure_path=sdf_full_path,smiles=smiles, fingerprint_type=FingerprintType.ECFP)
                ecfp_dict[inchi_key] = fingerprint_ecfp
                fingerprint_maccs = get_fingerprint(aDesc,mol_structure_path=sdf_full_path,smiles=smiles, fingerprint_type=FingerprintType.MACCSFP)
                maccsfp_dict[inchi_key] = fingerprint_maccs
                fingerprint_pfp = get_fingerprint(aDesc,mol_structure_path=sdf_full_path,smiles=smiles, fingerprint_type=FingerprintType.PFP)
                pfp_dict[inchi_key] = fingerprint_pfp
                try:
                    fingerprint_morgan = get_morgan_fingerprint_rdkit(chemicalStructureFile=sdf_full_path,smiles=smiles)
                except Exception as e: 
                    fingerprint_morgan = "NA"
                morganfp_dict[inchi_key] = fingerprint_morgan
                
                vector_fingerprints = generate_vector_fingerprints(aDesc,mol_structure_path=sdf_full_path,smiles=smiles)
                vector_fingerprints_dict[inchi_key] = vector_fingerprints
            
            
            partialDictDescriptorsAndFingerprintsRow['ECFP'] = fingerprint_ecfp
            partialDictDescriptorsAndFingerprintsRow['MACCSFP'] = fingerprint_maccs
            partialDictDescriptorsAndFingerprintsRow['PFP'] = fingerprint_pfp
            
            partialDictDescriptorsAndFingerprintsRow['MorganFP'] = fingerprint_morgan
            writerDescriptorsAndFingerprints.writerow(partialDictDescriptorsAndFingerprintsRow)
            
            partialDictFP = row.copy()
            for i in range(0,NUMBER_FPVALUES):
                header_name = "V" + str(i+1)
                partialDictFP[header_name] = vector_fingerprints[i]
                
            writerFingerprintsVectorized.writerow(partialDictFP)
            

