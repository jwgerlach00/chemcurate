import json
import duckdb
import pyarrow as pa
import copy
import zipfile
from typing import List, Dict, Iterable, Generator
import gzip
import os
from itertools import zip_longest
from pprint import pprint
from tqdm import tqdm
from rdkit.Chem import PandasTools
import naclo
import sqlite3
from ABCChemDB import ABCChemDB


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
        
        self.__possible_columns = [ # sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn
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
        
        ########## Connect to DB ##########
        self.conn = sqlite3.connect('pubchem_temp.db')
        self.cur = self.conn.cursor()
        
        ########## Create tables ##########
        self.cur.execute('DROP TABLE IF EXISTS substance') # remove substance table if it exists to avoid duplication
        self.build_substance_table_and_add_to_db()
        
        # self.read_whole_dir()
        # # print('o')
        # loader = self.read_zip_dir('/Users/collabpharma/Desktop/JSON/0434001_0435000.zip', verbose=True, batch_size=3)
        # for batch in loader:
        #         for json in batch:
        #             table = self.format_file(json)
        
    def build(self) -> None:
        ########## Build substance table ##########
        ########## Build bioassay table ##########
        ########## Join substances to bioassay table ##########
        ########## Remove substance table ##########
        return None # NOTE: Maybe return a connection to the DB?
    
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
        return copy.copy(self.__unit_map)

    @property
    def _activity_map(self) -> Dict[int, str]:
        """
        Sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn

        :return: Mapper from INTEGER code to STRING activity class
        :rtype: Dict[int, str]
        """
        return copy.copy(self.__activity_map)

    @property
    def _results_schema(self) -> pa.lib.Schema:
        return copy.copy(self.__results_schema)
    
    @property
    def _possible_columns(self) -> List[str]:
        """
        Sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn

        :return: List of names of possibe non-tid columns
        :rtype: List[str]
        """
        return copy.copy(self.__possible_columns)
    
    @staticmethod
    def _batch(iterable: Iterable[str], n: int = 1) -> Generator[List[str], None, None]:
        args = [iter(iterable)] * n
        yield from ([str(x) for x in tup if x is not None] for tup in zip_longest(*args))
    
    @staticmethod
    def read_zip_dir(zip_dir_path:str, verbose:bool=False, batch_size:int=10) -> List[dict]:
        """
        Reads all .json.gz files in a .zip directory and returns a list of dict objects

        :param zip_dir_path: Path to .zip directory containing .json.gz files
        :type zip_dir_path: str
        :param verbose: Whether to print filenames, defaults to False
        :type verbose: bool, optional
        :param batch_size: Number of .json.gz files to load into memory at once, defaults to 10
        :type batch_size: int, optional
        :yield: List of dict objects containing the loaded JSON data
        :rtype: Iterator[List[dict]]
        """
        # print(zip_dir_path)
        with zipfile.ZipFile(zip_dir_path, 'r') as zip_ref:
            # print(zip_ref)
            for batch in PubChemDB._batch(zip_ref.namelist(), n=batch_size):
                jsons = []
                # print(batch)
                # Iterate over .zip zipped directory
                for filename in batch:
                    if filename.endswith('.json.gz'):
                        # Read .json.gz zipped file
                        with zip_ref.open(filename, 'r') as file:
                            contents = gzip.decompress(file.read())
                            contents_str = contents.decode('utf-8')
                            # Load the JSON data into a Python object
                            if verbose: # Print filename
                                print(f'LOADING: {filename}')
                            jsons.append(json.loads(contents_str))
                yield jsons
    
    def read_whole_dir(self) -> None:
        """
        Iterates through entire directory OF ZIP DIRECTORIES, each containing .json.gz files and stores in database
        """
        for zip_dir in tqdm(os.listdir(self.json_dir_path)[:10]): # NOTE: Limit to just first for testing
            loader = self.read_zip_dir(os.path.join(self.json_dir_path, zip_dir), verbose=False, batch_size=3)
            
            for batch in loader:
                for json in batch:
                    table = self.format_file(json)
                    
                    # TODO: Add functionality to store in database
                    
    @staticmethod
    def _get_raw_data_from_file_json(file_json:dict) -> dict:
        """Gets deepcopy of data dict from JSON dict tree. Raises error if data is not present

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
        """Gets deepcopy of results dict from JSON dict tree. Raises error if results is not present

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
    def _format_raw_data_into_pa_table(raw_data:dict) -> pa.lib.Table:
        """_summary_

        :param data: _description_
        :type data: dict
        :return: _description_
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
        """Adds units to column names according to self.unit_map. ex: 'IC50' -> 'IC50, M'

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
        
    def format_file(self, file_json:dict) -> pa.lib.Table:
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
        except KeyError: # no data in file_json
            return None
        
        ########## Format raw data and results to PyArrow tables ##########
        data_table = PubChemDB._format_raw_data_into_pa_table(raw_data)
        results_table = pa.Table.from_pylist(raw_results, schema=self._results_schema)
        
        ########## Add units to each column name, ex: value --> value (unit) ##########
        results_table = PubChemDB._add_units_to_results_table_col_names(results_table)
            
        ########## Replace integer codes in data_table w/ names in results_table by joining on TIDs ##########
        non_tid_names = [i for i in self._possible_columns if i in data_table.column_names]
        tid_integer_codes = [int(x) for x in data_table.column_names if x not in non_tid_names]
        tid_names = [x[1] for x in
                     duckdb.sql(f'SELECT * FROM results_table WHERE tid IN {tuple(tid_integer_codes)}').fetchall()]
        data_table = data_table.rename_columns(non_tid_names + tid_names)

        ########## Convert integer coded activity to strings w/ meaning (same format as PubChem) ##########
        activity_col = pa.array(duckdb.sql(self.activity_map_sql_query).fetchall()).flatten()
        data_table = data_table.append_column('Activity', activity_col)
        # Remove old integer activity 'outcome' column
        data_table = data_table.remove_column(data_table.column_names.index('outcome'))

        return data_table
    
    def build_substance_table_and_add_to_db(self):
        for filename in os.listdir(self.sdf_dir_path):
            if filename.endswith('.sdf.gz'):
                df = PandasTools.LoadSDF(os.path.join(self.sdf_dir_path, filename))
                # Delete unnecessary columns
                df = df[['PUBCHEM_SUBSTANCE_ID', 'ROMol']]
                # Add SMILES from MOLs
                df = naclo.dataframes.df_mols_2_smiles(df, 'ROMol', 'SMILES')
                # Remove MOLs, only keep SMILES
                df.drop('ROMol', axis=1, inplace=True)
                
                # TODO: Add df to DB
                df.to_sql('substance', self.conn, if_exists='append', index=False)
    
if __name__ == '__main__':
    PubChemDB('/Users/collabpharma/Desktop/JSON',
              '/Users/collabpharma/Desktop/SDF')#/0001001_0002000.zip'
