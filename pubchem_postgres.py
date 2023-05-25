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
import psycopg2
from functools import wraps
from rdkit import RDLogger
from rdkit.rdBase import LogStatus as RDLogStatus
from sqlalchemy import create_engine


def rdkit_stfu(func):
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



class PubChemDB(ABCChemDB):
    def __init__(self, bioassay_json_dir_path:str, substance_sdf_dir_path:str) -> None:
        self.json_dir_path = bioassay_json_dir_path
        self.sdf_dir_path = substance_sdf_dir_path
        
        host = 'localhost'
        database = 'pubchem'
        user = 'postgres'
        password = 'Coll@bor@tions2020'
        port = '5432'

        
        ########## Connect to DB ##########
        self.connection = psycopg2.connect(host=host, database=database, user=user, password=password)
        self.cursor = self.connection.cursor()
        
        # Create a SQLAlchemy engine object
        self.engine = create_engine(f'postgresql+psycopg2://{user}:{password}@{host}:{port}/{database}')
        
        # self.cursor.execute("SELECT table_name FROM information_schema.tables WHERE table_schema='public';")

        # Fetch all the results
        # table_names = self.cursor.fetchall()
        # print(table_names)
        self.cursor.execute('SELECT * FROM substance')
        print(self.cursor.fetchall())
        
    def build(self):
        self._repopulate_substance_table()
        pass
    
    @rdkit_stfu
    def _repopulate_substance_table(self) -> None:
        ########## Remove existing tables ##########
        self.cursor.execute('DELETE FROM substance') # remove substance table if it exists to avoid duplication
        self.cursor.execute('DELETE FROM substance_errors')

        for filename in os.listdir(self.sdf_dir_path):
            if filename.endswith('.sdf.gz'):
                print(filename)
                try:
                    df = PandasTools.LoadSDF(os.path.join(self.sdf_dir_path, filename))
                    # Delete unnecessary columns
                    df = df[['PUBCHEM_SUBSTANCE_ID', 'ROMol']]
                    # Add SMILES from MOLs
                    df = naclo.dataframes.df_mols_2_smiles(df, 'ROMol', 'smiles')
                    # Remove MOLs, only keep SMILES
                    df = df.drop('ROMol', axis=1)
                    
                    # Insert each row into table
                    for _, row in df.iterrows():
                        self.cursor.execute('INSERT INTO substance (substance_id, smiles) VALUES (%s, %s)',
                                            (row['PUBCHEM_SUBSTANCE_ID'], row['smiles']))
                    
                    # Commit all executions
                    self.connection.commit()
                except Exception as e:
                    self.cursor.execute(f'INSERT INTO substance_errors (filename, error_message) VALUES ({0}, {1})'.format(
                        filename.split('.')[0], str(e)
                    ))
                    self.connection.commit()
        


if __name__ == '__main__':
    pc_db = PubChemDB('/Users/collabpharma/Desktop/JSON', '/Users/collabpharma/Desktop/SDF')
    pc_db.build()