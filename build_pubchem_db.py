import json
import duckdb
import pyarrow as pa
import copy
import zipfile
import tempfile
from typing import List, Dict, Iterable, Generator
import gzip
import os
import pyarrow.compute as pc
from itertools import zip_longest


class PubChemDB:
    def __init__(self, path):
        self.dir_path = path
        # self.file_jsons = self.read_dir(path)
        
        self.__unit_map = {
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
        
        self.__activity_map = {
            1: 'Inactive',
            2: 'Active',
            3: 'Inconclusive'
        }
        
        self.__results_schema = pa.schema([
            pa.field('tid', pa.int64()),
            pa.field('name', pa.string()),
            pa.field('unit', pa.int64()),
            pa.field('sunit', pa.string())
        ])
        
        self.read_whole_dir()
    
    @property  
    def unit_map(self) -> Dict[int, str]:
        return copy.deepcopy(self.__unit_map)

    @property
    def activity_map(self) -> Dict[int, str]:
        return copy.deepcopy(self.__activity_map)

    @property
    def results_schema(self) -> pa.lib.Schema:
        return copy.deepcopy(self.__results_schema)
    
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
        with zipfile.ZipFile(zip_dir_path, 'r') as zip_ref:
            for batch in PubChemDB.batch(zip_ref.namelist(), n=batch_size):
                jsons = []
                print(batch)
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
        for zip_dir in os.listdir(self.dir_path)[:1]: # NOTE: Limit to just first for testing
            loader = self.read_zip_dir(os.path.join(self.dir_path, zip_dir), verbose=False, batch_size=3)
            
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
        # Extract from file data
        data_copy = copy.deepcopy(file_json['PC_AssaySubmit']['data'])
        
        ########## Format data_copy ##########
        pylist = []
        for sid_entry in data_copy:
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
        exclude_names = [i for i in ['sid', 'version', 'outcome', 'rank', 'comment'] if i in data_table.column_names] # non-TID column names
        # List of TID column names in data_table
        old_names = [int(x) for x in data_table.column_names if x not in exclude_names]

        # List of names mapped to TIDs in results_table
        print(self.results_schema)
        results_table = pa.Table.from_pylist(file_json['PC_AssaySubmit']['assay']['descr']['results'], schema=self.results_schema)
        
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

        new_names = [x[1] for x in duckdb.sql(f'SELECT * FROM results_table WHERE tid IN {tuple(old_names)}').fetchall()]
        
        data_table = data_table.rename_columns(exclude_names + new_names)

        ########## Convert 1, 2, 3 coded activity to strings w/ meaning (same format as PubChem) ##########
        # If outcome is not able to be cast to an integer, set to NULL. If it is not 1, 2, or 3 set to NULL
        query = f'''
        SELECT
        CASE
            WHEN TRY_CAST(outcome AS INTEGER) IS NULL THEN NULL
            WHEN outcome = 1 THEN '{self.activity_map[1]}'
            WHEN outcome = 2 THEN '{self.activity_map[2]}'
            WHEN outcome = 3 THEN '{self.activity_map[3]}'
            ELSE NULl
        END AS Activity
        FROM data_table;
        '''

        activity_col = pa.array(duckdb.sql(query).fetchall()).flatten()
        data_table = data_table.append_column('Activity', activity_col)
        # Remove old integer activity 'outcome' column
        data_table = data_table.remove_column(data_table.column_names.index('outcome'))
        return data_table
    
if __name__ == '__main__':
    PubChemDB('/Users/collabpharma/Desktop/JSON')#/0001001_0002000.zip')
    
