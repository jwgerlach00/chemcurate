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


class PubChemDB:
    def __init__(self, json_dir_path:str, sdf_dir_path:str):
        self.json_dir_path = json_dir_path
        self.sdf_dir_path = sdf_dir_path
        
        self.__unit_map = { # Sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn
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
        
        self.__activity_map = { # Sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn
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
                WHEN outcome = 1 THEN '{self.activity_map[1]}'
                WHEN outcome = 2 THEN '{self.activity_map[2]}'
                WHEN outcome = 3 THEN '{self.activity_map[3]}'
                WHEN outcome = 4 THEN '{self.activity_map[4]}'
                WHEN outcome = 5 THEN '{self.activity_map[5]}'
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
        
        self.__possible_columns = [ # Sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn
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
        
        self.conn = sqlite3.connect('pubchem_temp.db')
        self.cur = self.conn.cursor()
        
        self.cur.execute('DROP TABLE IF EXISTS substance')
        self.build_substance_db()
        
        # self.read_whole_dir()
        # # print('o')
        # loader = self.read_zip_dir('/Users/collabpharma/Desktop/JSON/0434001_0435000.zip', verbose=True, batch_size=3)
        # for batch in loader:
        #         for json in batch:
        #             table = self.format_file(json)
    
    @property  
    def unit_map(self) -> Dict[int, str]:
        """
        Sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn

        :return: Mapper from INTEGER code to STRING unit of measurement
        :rtype: Dict[int, str]
        """
        return copy.copy(self.__unit_map)

    @property
    def activity_map(self) -> Dict[int, str]:
        """
        Sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn

        :return: Mapper from INTEGER code to STRING activity class
        :rtype: Dict[int, str]
        """
        return copy.copy(self.__activity_map)

    @property
    def results_schema(self) -> pa.lib.Schema:
        return copy.copy(self.__results_schema)
    
    @property
    def possible_columns(self) -> List[str]:
        """
        Sourced from: https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/pcassay2.asn

        :return: List of names of possibe non-tid columns
        :rtype: List[str]
        """
        return copy.copy(self.__possible_columns)
    
    @staticmethod
    def batch(iterable: Iterable[str], n: int = 1) -> Generator[List[str], None, None]:
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
            for batch in PubChemDB.batch(zip_ref.namelist(), n=batch_size):
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
        
    def format_file(self, file_json:dict) -> pa.lib.Table:
        """
        Formats a single PubChem FTP JSON (unzipped) file into a pyarrow table similar to that found in PubChem Web

        :param file_json: Loaded JSON file
        :type file_json: dict
        :return: Formatted pyarrow table
        :rtype: pa.lib.Table
        """
        # NO DATA, return null
        if 'data' not in file_json['PC_AssaySubmit'].keys() or \
            'results' not in file_json['PC_AssaySubmit']['assay']['descr'].keys():
            return None

        # Extract from file data
        data_copy = copy.deepcopy(file_json['PC_AssaySubmit']['data'])
        
        ########## Format data_copy ##########
        pylist = []
        for sid_entry in data_copy:
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
            
        data_table = pa.Table.from_pylist(pylist) # convert to table
            
        ########## Join data_copy w/ results_table on TIDs ##########
        exclude_names = [i for i in self.possible_columns if i in data_table.column_names] # non-TID column names
        # List of TID column names in data_table
        old_names = [int(x) for x in data_table.column_names if x not in exclude_names]
        # List of names mapped to TIDs in results_table
        results_table = pa.Table.from_pylist(file_json['PC_AssaySubmit']['assay']['descr']['results'],
                                             schema=self.results_schema)
        
        ########## Add units to name in results_table ##########
        results_dict = results_table.to_pydict()

        for i in range(len(results_table)):
            unit = results_dict['unit'][i]
            sunit = results_dict['sunit'][i]
            name = results_dict['name'][i]

            if unit is not None:
                results_dict['name'][i] = f'{name}, {self.unit_map[unit]}'
            elif sunit is not None:
                results_dict['name'][i] = f'{name}, {sunit}'

        results_table = pa.Table.from_pydict(results_dict)

        new_names = [x[1] for x in
                     duckdb.sql(f'SELECT * FROM results_table WHERE tid IN {tuple(old_names)}').fetchall()]
        
        data_table = data_table.rename_columns(exclude_names + new_names)

        ########## Convert integer coded activity to strings w/ meaning (same format as PubChem) ##########
        activity_col = pa.array(duckdb.sql(self.activity_map_sql_query).fetchall()).flatten()
        data_table = data_table.append_column('Activity', activity_col)
        # Remove old integer activity 'outcome' column
        data_table = data_table.remove_column(data_table.column_names.index('outcome'))
        print(data_table.column_names)
        return data_table
    
    def build_substance_db(self):
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