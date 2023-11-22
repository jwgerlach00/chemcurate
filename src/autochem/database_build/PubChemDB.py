import json
import pandas as pd
import copy
import zipfile
from typing import List, Dict, Optional, Generator
import gzip
import os
from tqdm import tqdm
from rdkit.Chem import PandasTools
import naclo
from __ABCChemDB import __ABCChemDB
import psycopg2
from functools import wraps
from rdkit import RDLogger
from rdkit.rdBase import LogStatus as RDLogStatus
import numpy as np
import json
import requests
import time
import logging


def _rdkit_stfu(func):
    """ Decorator function to silence rdkit warnings during evaluation of a function. """
    @wraps(func)
    def wrapper(*args, **kwargs):
        if RDLogStatus().count('disabled') < 4:
            RDLogger.DisableLog('rdApp.*') # rdkit is too verbose; this shuts it up
            returned_stuff = func(*args, **kwargs)
            RDLogger.EnableLog('rdApp.*') # rip off the duck tape; rdkit can speak again!
        else:
            # If rdkit log is already disabled, do not reenable it
            returned_stuff = func(*args, **kwargs)
        return returned_stuff
    return wrapper


class PubChemDB(__ABCChemDB):
    def __init__(self, bioassay_json_dir_path:str, substance_sdf_dir_path:str, protein2xrefs_path:str) -> None:
        self.bioassay_json_dir_path = bioassay_json_dir_path
        self.substance_sdf_dir_path = substance_sdf_dir_path
        self.protein2xrefs_path = protein2xrefs_path
        
        self.__activity_outcome_map = { # sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn
            1: 'inactive',
            2: 'active',
            3: 'inconclusive',
            4: 'unspecified',
            5: 'probe'
        }
        
        self.__possible_non_tid_bioassay_columns = [ # sourced from: \
            # https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn
            'sid',
            'sid_source', # source says "sid-source" but it is actually "sid_source"
            'version',
            'comment',
            'outcome',
            'rank',
            'data',
            'url',
            'xref',
            'date'
        ]
        
        self.__unit_map = { # sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn
            1: 'ppt',
            2: 'ppm',
            3: 'ppb',
            4: 'mm',
            5: 'um',
            6: 'nm',
            7: 'pm',
            8: 'fm',
            9: 'mgml',
            10: 'ugml',
            11: 'ngml',
            12: 'pgml',
            13: 'fgml',
            14: 'm',
            15: 'percent',
            16: 'ratio',
            17: 'sec',
            18: 'rsec',
            19: 'min',
            20: 'rmin',
            21: 'day',
            22: 'rday',
            23: 'ml-min-kg',
            24: 'l-kg',
            25: 'hr-ng-ml',
            26: 'cm-sec',
            27: 'mg-kg',
            254: 'none',
            255: 'unspecified'
        }
        
        ########## Connect to DB ##########
        host = 'localhost'
        database = 'pubchem'
        user = 'postgres'
        password = 'Coll@bor@tions2020'
        # port = '5432'
        self.connection = psycopg2.connect(host=host, database=database, user=user, password=password)
        self.cursor = self.connection.cursor()
        
    @property
    def _activity_outcome_map(self) -> Dict[int, str]:
        """
        Sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn

        :return: Mapper from INTEGER code to STRING activity class
        :rtype: Dict[int, str]d
        """
        return copy.deepcopy(self.__activity_outcome_map)
    
    @property
    def _possible_non_tid_bioassay_columns(self) -> List[str]:
        """
        Sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn

        :return: List of names of possibe non-tid columns
        :rtype: List[str]
        """
        return copy.deepcopy(self.__possible_non_tid_bioassay_columns)
    
    @property  
    def _unit_map(self) -> Dict[int, str]:
        """
        Sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn

        :return: Mapper from INTEGER code to STRING unit of measurement
        :rtype: Dict[int, str]
        """
        return copy.deepcopy(self.__unit_map)
    
    def uniprot_id_assay_id_map_relation(self):
        logging.basicConfig(filename='error_log.txt', level=logging.ERROR)
        
        print('Rebuilding uniprot_id_assay_id_map relation...')
        print('...clearing existing data...')
        try:
            self.cursor.execute('DROP TABLE BioassayToUniprot')
            self.cursor.execute('CREATE TABLE BioassayToUniprot (bioassayID INT, uniprotID VARCHAR(20), PRIMARY KEY \
                (bioassayID, uniprotID))')
            self.connection.commit()
        except Exception as e:
            print(f'ERROR: {e}, rolling back and exiting...')
            self.connection.rollback()
            exit()
        print('...getting all uniprot_ids...')
        self.cursor.execute('SELECT uniprot_id FROM target WHERE uniprot_id IS NOT NULL')
        uniprot_ids = self.cursor.fetchall()
        
        print('...getting all assay_ids')
        for uniprot_id in uniprot_ids: # probably a problem with the single transaction lasting too long
            t_0 = time.time()
            uniprot_id = uniprot_id[0]
            incomplete_flag = True
            while incomplete_flag:
                try:
                    _url_stem = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
                    url = f'{_url_stem}/bioassay/target/ProteinName/{uniprot_id}/aids/JSON'
                    
                    json_response = requests.get(url).json()
                    if 'IdentifierList' not in json_response:
                        continue
                    
                    assay_ids = requests.get(url).json()['IdentifierList']['AID']
                    
                    for assay_id in assay_ids:
                        try:
                            self.cursor.execute('INSERT INTO BioassayToUniprot VALUES (%s, %s)', (assay_id, uniprot_id))
                        except psycopg2.IntegrityError as e:
                            # Check if the error is a duplicate key violation
                            if 'duplicate key value violates unique constraint' in str(e):
                                print(f'Skipping duplicate row: bioassay_id={assay_id}, uniprot_id={uniprot_id}')
                            else:
                                # For other IntegrityError cases, print the error and rollback
                                print(f'Rolling back due to error: {e}')
                                self.connection.rollback()
                        else:
                            # Commit only if there was no exception during execution
                            self.connection.commit()
                except Exception as e:
                    print(f'{uniprot_id}, ERROR: {e}')
                    logging.error(f"{uniprot_id}, {str(e)}")
                    self.connection.rollback()
                else:
                    logging.error(f"{uniprot_id}")
                    incomplete_flag = False
                
            t_diff = time.time() - t_0
            if t_diff < 0.20:
                time.sleep(0.20 - t_diff)  # NOTE: Required because PubChem API has a limit of 5 requests per second
        
    def build(self):
        self.repopulate_protein_target_table() # NOTE: done
        self.repopulate_substance_table()
        self.repopulate_bioassay_table()
        pass
    
    def repopulate_protein_target_table(self):
        # Clear existing data to avoid duplication
        self.cursor.execute('DELETE FROM target')
        
        df = pd.read_csv(self.protein2xrefs_path, delimiter='\t')
        
        # Replace w/ None instead of NaN because NaN is a string while None is recognized as NULL
        df = df.replace({np.nan: None})
        
        # There are some duplicate protein accessions w/ different GeneIDs, this won't work as protein_accession is \
            # the primary key
        df = df.drop_duplicates(subset=['ProteinAccession'], keep='first')
        
        # Insert each row into table
        for _, row in df.iterrows():
            self.cursor.execute('INSERT INTO target (protein_accession, uniprot_id) VALUES (%s, %s)',
                                (row['ProteinAccession'], row['UniProt']))

        self.connection.commit()
        
    @_rdkit_stfu
    def repopulate_substance_table(self) -> None:
        # Clear existing data to avoid duplication
        self.cursor.execute('DELETE FROM substance')
        self.cursor.execute('DELETE FROM substance_errors')
        self.connection.commit()

        for filename in tqdm(os.listdir(self.substance_sdf_dir_path)):
            if filename.endswith('.sdf.gz'):
                try:
                    df = PandasTools.LoadSDF(os.path.join(self.substance_sdf_dir_path, filename))
                    # Delete unnecessary columns
                    df = df[['PUBCHEM_SUBSTANCE_ID', 'ROMol']]
                    # Add SMILES from MOLs
                    df = naclo.dataframes.df_mols_2_smiles(df, 'ROMol', 'smiles')
                    # Remove MOLs, only keep SMILES
                    df = df.drop('ROMol', axis=1)
                    
                    # Replace w/ None instead of NaN because NaN is a string while None is recognized as NULL
                    df = df.replace({np.nan: None})
                    
                    # Insert each row into table
                    for _, row in df.iterrows():
                        self.cursor.execute('INSERT INTO substance (substance_id, smiles) VALUES (%s, %s)',
                                            (row['PUBCHEM_SUBSTANCE_ID'], row['smiles']))

                except Exception as e:
                    self.cursor.execute('INSERT INTO substance_errors (filename, error_message) VALUES (%s, %s)',
                                        (filename.split('.')[0], str(e)))
                
                # Commit executions
                self.connection.commit()
                
    def repopulate_bioassay_table(self, protein_only:bool=True) -> None:
        ########## Clear existing data to avoid duplication ##########
        self.cursor.execute('DELETE FROM bioassay')
        self.connection.commit()
        
        for zip_dir in tqdm([file for file in os.listdir(self.bioassay_json_dir_path) if not 
                             file.startswith('README')]): # NOTE: Index os.listdir to limit number of files for testing
            
            loader = PubChemDB._bioassay_zip_dir_loader(os.path.join(self.bioassay_json_dir_path, zip_dir),
                                                        print_filename=False)
            
            # Skip if loader is None meaning the zip_dir failed to load
            if loader == None:
                continue
            
            for bioassay_json in loader:
                # Skip if from chembl
                if bioassay_json['PC_AssaySubmit']['assay']['descr']['aid_source']['db']['name'].lower() == 'chembl':
                    continue
                
                self._protein_only_add_entry_to_bioassay_table(bioassay_json)
    
    @staticmethod
    def _bioassay_zip_dir_loader(zip_dir_path:str, print_filename:bool=False) -> Optional[Generator]:
        """
        Sequentially reads .json.gz files in a .zip directory into dict object. Yields one at a time.
        
        NOTE: This will read from a zip directory with .json.gz children, data comes directly from FTP as a zip
        directory containing ZIP DIRECTORIES, each of which contain .json.gz files.

        :param zip_dir_path: Path to .zip directory containing .json.gz files
        :type zip_dir_path: str
        :param print_filename: Whether to print filenames, defaults to False
        :type print_filename: bool, optional
        :yield: Dict object containing loaded JSON data
        :rtype: Iterator[List[dict]]
        """
        try:
            with zipfile.ZipFile(zip_dir_path, 'r') as zip_ref:
                # Iterate over .zip zipped directory
                filenames = [filename for filename in zip_ref.namelist() if filename.endswith('.json.gz')]

                for filename in filenames:
                    # Read .json.gz zipped file
                    with zip_ref.open(filename, 'r') as file:
                        contents = gzip.decompress(file.read())
                        try:
                            contents_str = contents.decode('utf-8')
                        except Exception as e:
                            continue
                        # Load the JSON data into a Python object
                        if print_filename: # Print filename
                            print(f'LOADING: {filename}')
        
                    yield json.loads(contents_str)
        except zipfile.BadZipFile as e:
            print(f'Failed to load zipdir: {zip_dir_path}') # likely means the zipdir is corrupted from FTP download
            return None
                
    def _protein_only_add_entry_to_bioassay_table(self, json_bioassay:dict) -> None:
        '''TODO: fails if not protein only'''
        bioassay_id = json_bioassay['PC_AssaySubmit']['assay']['descr']['aid']['id']
        
        ########## Determine protein accession from bioassay target info ##########
        target_info = PubChemDB._get_bioassay_target_info_if_exists(json_bioassay)
        
        if target_info and 'protein_accession' in target_info['mol_id']: # NOTE: target_info could be NoneType or \
            # empty object
            protein_accession = target_info['mol_id']['protein_accession']
            protein_accession = protein_accession.split('.')[0] # remove version number if present
        else:
            protein_accession = None
        
        ########## Format bioassay data ########## 
        if 'data' in json_bioassay['PC_AssaySubmit']:
            formatted_bioassay_data = self._reformat_bioassay_data(json_bioassay)
        else: # no data for the bioassay
            formatted_bioassay_data = None
        
        ########## Store in database ##########
        try:
            self.cursor.execute('INSERT INTO bioassay (bioassay_id, protein_accession, assay_data) VALUES (%s, %s, %s)',
                                (bioassay_id, protein_accession, json.dumps(formatted_bioassay_data)))
        except psycopg2.errors.ForeignKeyViolation: # the protein accession is not in the protein table
            self.connection.rollback() # rollback to previous commit due to foreign key violation
            self.cursor.execute('INSERT INTO bioassay (bioassay_id) VALUES (%s)', (bioassay_id,))

        self.connection.commit()
        
    @staticmethod
    def _get_bioassay_target_info_if_exists(json_bioassay:dict) -> Optional[dict]:
        try:
            target_info = json_bioassay['PC_AssaySubmit']['assay']['descr']['target'][0] # NOTE: it comes as a list of \
                # one element for some reason
            return target_info
        except (KeyError, IndexError): # KeyError for ['target'] or IndexError for [0], assume no target data
            return None
    
    def _reformat_bioassay_data(self, json_bioassay:dict) -> dict:
        if 'results' in json_bioassay['PC_AssaySubmit']['assay']['descr']:
            tid_to_activity_name_map = {}
            for item in json_bioassay['PC_AssaySubmit']['assay']['descr']['results']:
                if 'unit' in item:
                    tid_to_activity_name_map[item['tid']] = f"{item['name']} ({self._unit_map[item['unit']]})"
                else:
                    tid_to_activity_name_map[item['tid']] = '' # no unit --> this will show up as a blank: "()"
        else:
            return [] # no TIDs therefore assume no assay data for the bioassay
        
        reformatted_data = []
        for sid_entry in json_bioassay['PC_AssaySubmit']['data']:
            if 'data' in sid_entry: # clean up data, else there is no activity data but the structure is still added \
                # to reformatted_data
                sid_results = sid_entry.pop('data')
                tid_data = {tid_entry['tid']: list(tid_entry['value'].values())[0] for tid_entry in sid_results} # \
                    # list(tid_entry['value'].values())[0] converts a dictionary with various keys such as 'sval' or \
                        # 'fval' to a single value
                sid_entry.update(tid_data)
            
            # Decode activity outcome from integers to strings
            sid_entry['outcome'] = self._activity_outcome_map[sid_entry['outcome']]
            
            # Decode TID to activity name w/ units
            sid_entry = {(tid_to_activity_name_map[key] if key not in self._possible_non_tid_bioassay_columns else key):
                value for key, value in sid_entry.items()}

            reformatted_data.append(sid_entry)
        return reformatted_data

if __name__ == '__main__':
    pc_db = PubChemDB(
        '/Users/collabpharma/Desktop/JSON',
        '/Users/collabpharma/Desktop/SDF',
        '/Users/collabpharma/Desktop/protein2xrefs'
    )
    pc_db.uniprot_id_assay_id_map_relation()