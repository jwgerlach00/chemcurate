import json
import duckdb
import pandas as pd
import pyarrow as pa
import copy
import zipfile
from typing import List, Dict
import gzip
import os
from tqdm import tqdm
from rdkit.Chem import PandasTools
import naclo
import sqlite3
from ABCChemDB import ABCChemDB
import re
from pprint import pprint


def sanitize_column_name(column_name:str, general_replacement:str='_'):
    # Replace < with 'less_than'
    sanitized_name = re.sub(r'<', 'less_than', column_name)
    
    # Replace > with 'greater_than'
    sanitized_name = re.sub(r'>', 'greater_than', sanitized_name)
    
    # Replace = with 'equals'
    sanitized_name = re.sub(r'=', 'equals', sanitized_name)
    
    # Replace other invalid characters with the specified replacement text
    sanitized_name = re.sub(r'[^a-zA-Z0-9_, ]', general_replacement, sanitized_name)
    
    # Ensure the column name starts with a letter
    if not sanitized_name[0].isalpha():
        sanitized_name = general_replacement + sanitized_name
    
    return sanitized_name


def print_table(cursor, table_name:str, n_rows:int=100) -> None:
    cursor.execute(f'SELECT * FROM {table_name}')
    rows = cursor.fetchall()
    # Print the table contents
    for row in rows:
        print(row)


class PubChemDB(ABCChemDB):
    def __init__(self, bioassay_json_dir_path:str, substance_sdf_dir_path:str) -> None:
        self.json_dir_path = bioassay_json_dir_path
        self.sdf_dir_path = substance_sdf_dir_path
        
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
        
        self.__activity_map = { # sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn
            1: 'inactive',
            2: 'active',
            3: 'inconclusive',
            4: 'unspecified',
            5: 'probe'
        }
        
        # If outcome is not able to be cast to an integer, set to NULL. If it is not a registered INT set to NULL
        self.activity_map_sql_query = f'''
            SELECT
            CASE
                WHEN TRY_CAST(outcome AS INTEGER) IS NULL THEN NULL
                WHEN outcome = 1 THEN '{self._activity_map[1]}'
                WHEN outcome = 2 THEN '{self._activity_map[2]}'
                WHEN outcome = 3 THEN '{self._activity_map[3]}'
                WHEN outcome = 4 THEN '{self._activity_map[4]}'
                WHEN outcome = 5 THEN '{self._activity_map[5]}'
                ELSE NULl
            END AS Activity
            FROM data_table;
        '''
        
        self.__results_schema = pa.schema([
            pa.field('tid', pa.int64()),
            pa.field('name', pa.string()),
            pa.field('unit', pa.int64()),
            pa.field('sunit', pa.string())
        ])
        
        self.__possible_non_tid_bioassay_columns = [ # sourced from: \
            # https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn
            'sid',
            'sid-source',
            'version',
            'comment',
            'outcome',
            'rank',
            'data',
            'url',
            'xref',
            'date'
        ]
        
        self.__possible_target_keys = [ # sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn
            'name',
            'mol_id',
            'organism',
            'descr',
            'comment'
        ]
        
        ########## Connect to DB ##########
        self.conn = sqlite3.connect('pubchem_temp.db')
        self.cur = self.conn.cursor()
        # TODO: Functionality for building a new DB called 'pubchem_temp' or wtvr
    
    @property
    def connection(self) -> sqlite3.Connection:
        return self.conn
    
    @property  
    def _unit_map(self) -> Dict[int, str]:
        """
        Sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn

        :return: Mapper from INTEGER code to STRING unit of measurement
        :rtype: Dict[int, str]
        """
        return copy.deepcopy(self.__unit_map)

    @property
    def _activity_map(self) -> Dict[int, str]:
        """
        Sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn

        :return: Mapper from INTEGER code to STRING activity class
        :rtype: Dict[int, str]
        """
        return copy.deepcopy(self.__activity_map)

    @property
    def _results_schema(self) -> pa.lib.Schema:
        return copy.deepcopy(self.__results_schema)
    
    @property
    def _possible_non_tid_bioassay_columns(self) -> List[str]:
        """
        Sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn

        :return: List of names of possibe non-tid columns
        :rtype: List[str]
        """
        return copy.deepcopy(self.__possible_non_tid_bioassay_columns)

    @property
    def _possible_target_keys(self) -> List[str]:
        """
        Sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn

        :return: List of keys that could be present in the raw target data
        :rtype: List[str]
        """
        return copy.deepcopy(self.__possible_target_keys)
        
    def build(self) -> None:
        ########## Build target tables ##########
        # self._build_protein_table()
        
        ########## Build substance table ##########
        # self._build_substance_table()
        
        ########## Build bioassay table ##########
        self._build_bioassay_table()
        
        # print_table(self.cur, 'bioassay')
        ########## Join substances to bioassay table ##########
        ########## Remove substance table ##########
        return None # NOTE: Maybe return a connection to the DB?

    def _build_substance_table(self) -> None:
        ########## Remove existing tables ##########
        self.cur.execute('DROP TABLE IF EXISTS substance') # remove substance table if it exists to avoid duplication
        self.cur.execute('DROP TABLE IF EXISTS substance_errors')
        ########## Create substance_errors table (NOTE: substance_table is created in \
            # build_substance_table_and_add_to_db()) ##########
        self.cur.execute('CREATE TABLE substance_errors (filename TEXT PRIMARY KEY, error TEXT)')
        self._build_substance_table_and_add_to_db()
        
    def _build_protein_table(self) -> None:
        path = '/Users/collabpharma/Desktop/protein2xrefs'
        df = pd.read_csv(path, delimiter='\t')
        df.to_sql('protein', self.conn, if_exists='replace', index=False)
        
    def _build_bioassay_table(self) -> None:
        self._clear_bioassay_target_relationship_table_and_specify_schema() # create table but don't yet populate
        
        # Remove existing tables
        self.cur.execute('DROP TABLE IF EXISTS bioassay')
        
        first_bioassay_flag_for_alter_table = True # if True, don't need to alter table to add new columns
            # initial schema = schema for first bioassay
        for zip_dir in tqdm(os.listdir(self.json_dir_path)): # NOTE: Limit to just first few for testing
            loader = PubChemDB._read_one_bioassay_zip_dir(os.path.join(self.json_dir_path, zip_dir),
                                                          print_filename=False)
            
            for json in loader: # (batch_size) bioassays at a time
                # for json in batch:
                table = self._format_bioassay_and_store_in_table(json)
                
                if not table: # table is None
                    continue
                
                if not first_bioassay_flag_for_alter_table:
                    existing_columns = [x[1] for x in self.conn.execute(f'PRAGMA table_info(bioassay)').fetchall()]
                    print(len(existing_columns))
                    
                    for column in table.column_names:
                        if column not in existing_columns:
                            column_type = table.schema.field(column).type
                            
                            # print(column)
                            # print(sanitize_column_name(column))
                            # print(table.to_pandas())
                            # print(f"ALTER TABLE bioassay ADD COLUMN '{column}' \
                            #     {column_type}")
                            self.conn.execute(f'ALTER TABLE bioassay ADD COLUMN "{column}" \
                                {column_type}')

                table.to_pandas().to_sql('bioassay', self.conn, if_exists='append', index=False)
                first_bioassay_flag_for_alter_table = False
                
                # self.cur.execute('SELECT * FROM bioassay')
                # rows = self.cur.fetchall()

                # # Print the table contents
                # for row in rows[:100]:
                #     print(row)   
                
                # TODO: Add targets as a column
                
                # TODO: Add functionality to store in database
                
    def _clear_bioassay_target_relationship_table_and_specify_schema(self):
        # Remove existing tables
        self.cur.execute('DROP TABLE IF EXISTS bioassay_target_relationship')
        create_table_sql = '''
            CREATE TABLE bioassay_target_relationship (
                bioassay_id INTEGER,
                uniprot_id TEXT
            )
        '''
        self.cur.execute(create_table_sql)
        self.conn.commit()
    
    def _build_substance_table_and_add_to_db(self):
        for filename in os.listdir(self.sdf_dir_path):
            if filename.endswith('.sdf.gz'):
                try:
                    df = PandasTools.LoadSDF(os.path.join(self.sdf_dir_path, filename))
                    # Delete unnecessary columns
                    df = df[['substance_id', 'ROMol']]
                    # Add SMILES from MOLs
                    df = naclo.dataframes.df_mols_2_smiles(df, 'ROMol', 'SMILES')
                    # Remove MOLs, only keep SMILES
                    df.drop('ROMol', axis=1, inplace=True)
                    
                    # TODO: Add df to DB
                    df.to_sql('substance', self.conn, if_exists='append', index=False)
                except Exception as e:
                    self.cur.execute('INSERT INTO substance_errors (filename, error) VALUES (?, ?)',
                                     (filename.split('.')[0], str(e)))
    
    @staticmethod
    def _read_one_bioassay_zip_dir(zip_dir_path:str, print_filename:bool=False) -> dict:
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
        with zipfile.ZipFile(zip_dir_path, 'r') as zip_ref:
            # Iterate over .zip zipped directory
            filenames = [filename for filename in zip_ref.namelist() if filename.endswith('.json.gz')]

            for filename in filenames:
                # Read .json.gz zipped file
                with zip_ref.open(filename, 'r') as file:
                    contents = gzip.decompress(file.read())
                    contents_str = contents.decode('utf-8')
                    # Load the JSON data into a Python object
                    if print_filename: # Print filename
                        print(f'LOADING: {filename}')
    
                yield json.loads(contents_str)
                
    def _format_bioassay_and_store_in_table(self, file_json:dict, protein_only:bool=True) -> pa.lib.Table:
        """
        Formats a single PubChem FTP JSON (unzipped) file into a pyarrow table similar to that found in PubChem Web

        :param file_json: Loaded JSON file
        :type file_json: dict
        :return: Formatted pyarrow table
        :rtype: pa.lib.Table
        """
        ########## Get raw data and results nested within file_json ##########
        try:
            raw_data = PubChemDB._get_raw_data_from_file_json(file_json)
            raw_results = PubChemDB._get_raw_results_from_file_json(file_json) # contains metadata
            formatted_target = self._format_raw_target(PubChemDB._get_raw_target_from_file_json(file_json))
        except KeyError: # no data in file_json
            return None
        
        ########## Check if target is protein ##########
        if protein_only:
            if 'target_protein_accession' not in formatted_target.keys(): # not protein
                return None
            else:
                bioassay_id = file_json['PC_AssaySubmit']['assay']['descr']['aid']['id']
                self._populate_bioassay_target_relationship_table(bioassay_id,
                                                                  formatted_target['target_protein_accession'])
                
        
        ########## Format raw data and results to PyArrow tables ##########
        data_table = PubChemDB._format_raw_data_into_pa_table(raw_data)
        results_table = pa.Table.from_pylist(raw_results, schema=self._results_schema)
        
        ########## Add units to each column name, ex: value --> value (unit) ##########
        results_table = self._add_units_to_results_table_col_names(results_table)
            
        ########## Replace integer codes in data_table w/ names in results_table by joining on TIDs ##########
        non_tid_names = [i for i in self._possible_non_tid_bioassay_columns if i in data_table.column_names]
        tid_integer_codes = [int(x) for x in data_table.column_names if x not in non_tid_names]
        tid_names = [x[1] for x in
                    duckdb.sql(f'SELECT * FROM results_table WHERE tid IN {tuple(tid_integer_codes)}').fetchall()]
        data_table = data_table.rename_columns(non_tid_names + tid_names)

        ########## Convert integer coded activity to strings w/ meaning (same format as PubChem) ##########
        activity_col = pa.array(duckdb.sql(self.activity_map_sql_query).fetchall()).flatten()
        data_table = data_table.append_column('Activity', activity_col)
        # Remove old integer activity 'outcome' column
        data_table = data_table.remove_column(data_table.column_names.index('outcome'))
        
        ########## Format target info and add as columns to data_table ##########
        # for key, value in self._format_raw_target(raw_target).items():
        #     column = pa.array([value] * len(data_table))
        #     data_table = data_table.add_column(data_table.num_columns, key, column)
        
        ########## Make the column names lowercase ##########
        columns_lower = [column.lower() for column in data_table.column_names]
        data_table = pa.Table.from_arrays(data_table.columns, schema=pa.schema([(column_lower, field.type) for
                                                                                column_lower, field in
                                                                                zip(columns_lower, data_table.schema)]))

        return data_table
    
    def _populate_bioassay_target_relationship_table(self, bioassay_id:int, uniprot_id:str) -> None:
        self.cur.execute('INSERT INTO bioassay_target_relationship (bioassay_id, uniprot_id) VALUES (?, ?)',
                         (bioassay_id, uniprot_id))
                    
    @staticmethod
    def _get_raw_data_from_file_json(file_json:dict) -> dict:
        """
        Gets deepcopy of data dict from JSON dict tree. Raises error if data is not present

        :param file_json: JSON file loaded to memory
        :type file_json: dict
        :raises KeyError: No data exists in file_json
        :return: Deepcopy of raw data dict from nested JSON
        :rtype: dict
        """
        if 'data' in file_json['PC_AssaySubmit'].keys():
            return copy.deepcopy(file_json['PC_AssaySubmit']['data'])
        else:
            raise KeyError('No data in file_json')
        
    @staticmethod
    def _get_raw_results_from_file_json(file_json:dict) -> dict:
        """
        Gets deepcopy of results dict from JSON dict tree. Raises error if results is not present

        :param file_json: JSON File loaded to memory
        :type file_json: dict
        :raises KeyError: No results exist in file_json
        :return: Deepcopy of raw results dict from nested JSON
        :rtype: dict
        """
        if 'results' in file_json['PC_AssaySubmit']['assay']['descr'].keys():
            return copy.deepcopy(file_json['PC_AssaySubmit']['assay']['descr']['results'])
        else:
            raise KeyError('No results in file_json')
    
    @staticmethod
    def _get_raw_target_from_file_json(file_json:dict):
        """
        Gets deepcopy of target dict from JSON dict tree. Raises error if target is not present

        :param file_json: JSON File loaded to memory
        :type file_json: dict
        :raises KeyError: No target exists in file_json
        :return: Deepcopy of raw target dict from nested JSON
        :rtype: dict
        """
        if 'target' in file_json['PC_AssaySubmit']['assay']['descr'].keys():
            try:
                return file_json['PC_AssaySubmit']['assay']['descr']['target'][0] # NOTE: it comes as a list of one \
            except Exception:
                print(file_json['PC_AssaySubmit']['assay']['descr']['target'])
                raise IndexError('ape')
                
                # element for some reason
        else:
            raise KeyError('No target in file_json')
        
    @staticmethod 
    def _format_raw_data_into_pa_table(raw_data:dict) -> pa.lib.Table:
        """
        _summary_

        :param data: Raw data info from JSON
        :type data: dict
        :return: Formatted table from the raw data
        :rtype: pa.lib.Table
        """
        pylist = []
        for sid_entry in raw_data:
            # NO DATA, skip
            if 'data' not in sid_entry.keys():
                continue
            sid_results = sid_entry.pop('data') # List of dicts

            # Extract useful data for each TID
            for i, tid_data in enumerate(sid_results):
                tid_data.update({str(tid_data.pop('tid')): # use string TID as dict key
                                list(tid_data.pop('value').values())[0]}) # list of one element
                sid_results[i] = tid_data # over-write

            # Convert list of dictionaries to single dict
            sid_results = {k: v for d in sid_results for k, v in d.items()}

            sid_entry.update(sid_results)
            pylist.append(sid_entry)
            
        return pa.Table.from_pylist(pylist) # convert to table
    
    def _add_units_to_results_table_col_names(self, results_table:pa.lib.Table) -> pa.lib.Table:
        """
        Adds units to column names according to self.unit_map. ex: 'IC50' -> 'IC50, M'

        :param results_table: PyArrow table containing PubChem results data
        :type results_table: pa.lib.Table
        :return: PyArrow table containing PubChem results data with units in column names
        :rtype: pa.lib.Table
        """
        results_dict = results_table.to_pydict()

        for i in range(len(results_table)):
            unit = results_dict['unit'][i]
            sunit = results_dict['sunit'][i]
            name = results_dict['name'][i]

            if unit is not None:
                results_dict['name'][i] = f'{name}, {self._unit_map[unit]}'
            elif sunit is not None:
                results_dict['name'][i] = f'{name}, {sunit}'

        return pa.Table.from_pydict(results_dict)
    
    def _format_raw_target(self, raw_target:dict) -> dict:
        """
        Un-nests the 'mol_id' data and renames some keys to make them appropriate column names

        :param raw_target: Raw target info from JSON, double-level dict
        :type raw_target: dict
        :return: Reorganized target info, single-level dict
        :rtype: dict
        """
        target_dict = {}
        
        for key in self._possible_target_keys:
            if key in raw_target.keys():
                if key == 'mol_id': # un-nest from double-level to single-level dict
                    for key, value in raw_target['mol_id'].items():
                        target_dict[f'target_{key}'] = value
                elif key == 'organism':
                    continue # skip, its annoying to store; TODO: figure out how to store this
                else:
                    target_dict[f'target_{key}'] = raw_target[key]
        
        return target_dict
        

if __name__ == '__main__':
    pc_db = PubChemDB('/Users/collabpharma/Desktop/JSON',
                      '/Users/collabpharma/Desktop/SDF')#/0001001_0002000.zip'
    pc_db.build()
