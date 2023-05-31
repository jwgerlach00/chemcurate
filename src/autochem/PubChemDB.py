from typing import List
import json
import gzip
import os
import tempfile
import zipfile
import shutil
import pandas as pd
from pprint import pprint
from autochem import PubChem


def read_json_gz(dir_path:str) -> List[dict]:
    data = []
    
    # Loop through each compressed folder in the directory and load the JSON data
    for foldername in os.listdir(dir_path):
        # if not foldername.endswith('.zip'):
        #     continue
    
        # Extract the compressed folder to a temporary directory
        tempdir = tempfile.mkdtemp()
        with zipfile.ZipFile(os.path.join(dir_path, foldername), 'r') as zip_ref:
            zip_ref.extractall(tempdir)
            
        # Loop through each JSON file in the temporary directory and load the data
        for filename in os.listdir(tempdir):
            if filename.endswith('.json.gz'):
                with gzip.open(os.path.join(tempdir, filename), 'rb') as file:
                    file_contents = file.read().decode('utf-8')
                    json_data = json.loads(file_contents)
                    data.append(json_data)
        # Delete the temporary directory
        shutil.rmtree(tempdir)
        
def read_one(dir_path:str) -> List[dict]:
    li = os.listdir(dir_path)
    
    tempdir = tempfile.mkdtemp()
    with zipfile.ZipFile(os.path.join(dir_path, '0001001_0002000.zip'), 'r') as zip_ref:
        for filename in zip_ref.namelist():
            if filename.endswith('.json.gz'):
                with zip_ref.open(filename, 'r') as file:
                    contents = gzip.decompress(file.read())
                    contents_str = contents.decode('utf-8')
                    # Load the JSON data into a Python object
                    print(filename)
                    data = json.loads(contents_str)
                    # pprint(data['PC_AssaySubmit'])
                    break
    #         zip_ref.extractall(tempdir)
    df = PubChem.get_concise_assay_dfs(['1001'])
    print(df.columns)
    
    # for filename in os.listdir(tempdir):
    #     print(filename)


if __name__ == '__main__':
    import sqlite3
    
    conn = sqlite3.connect('pubchem_temp.db')
    cursor = conn.cursor()
    
    path = '../../../Desktop/JSON'
    
    read_one(path)
    