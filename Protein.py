import requests
import pandas as pd
from typing import Union
import time


# def get_full_bioassay_record(assay_id:Union[int, str]):
#     url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/bioassay/aid/{assay_id}/record/CSV'
#     response = requests.get(url)
#     return response.text

# def get_assay_summary(assay_id:Union[int, str]):
#     url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/bioassay/aid/{assay_id}/summary/JSON'
#     response = requests.get(url)
#     return response.text


class Protein:
    def __init__(self, uniprot_id:str) -> None:
        self.__smiles_col_name = 'SMILES'
        
        self.uniprot = uniprot_id
        
        self.concise_assay_dfs = []
        self.sid_records = []
        self.cids = []
        self.smiles = []
        
        self.assay_ids = self.__get_assay_ids()
        self.__set_concise_assay_dfs()
        self.__set_sid_records()
        self.__cids_from_sid_records()
        self.__set_smiles()
        self.__add_smiles_to_dfs()
        
        self.df = pd.concat(self.concise_assay_dfs)
        self.df.reset_index(drop=True, inplace=True)
        
        smiles_col = self.df.pop(self.__smiles_col_name)
        self.df.insert(0, self.__smiles_col_name, smiles_col)
    
    @staticmethod
    def get_assay_ids(uniprot:str) -> list:
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/bioassay/target/ProteinName/{uniprot}/aids/JSON'
        response = requests.get(url)
        return response.json()['IdentifierList']['AID']
    
    def __get_assay_ids(self) -> list:
        return Protein.get_assay_ids(self.uniprot)
    
    @staticmethod
    def get_concise_assay_dfs(assay_id:Union[int, str]) -> pd.DataFrame:
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/bioassay/aid/{assay_id}/concise/JSON'
        response = requests.get(url)
        json = response.json()['Table']
        rows = [row['Cell'] for row in json['Row']]
        return pd.DataFrame(rows, columns=json['Columns']['Column'])

    def __set_concise_assay_dfs(self):
        for aid in self.assay_ids:
            self.concise_assay_dfs.append(Protein.get_concise_assay_dfs(aid))
            time.sleep(0.5)
            
    @staticmethod
    def get_sid_records(sids):
        sids_str = ','.join(sids)
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sids_str}/record/JSON'
        response = requests.get(url)
        return response.json()['PC_Substances']
    
    def __set_sid_records(self):
        for df in self.concise_assay_dfs:
            self.sid_records.append(Protein.get_sid_records(df['SID'].tolist()))
            time.sleep(0.5)
    
    @staticmethod
    def cids_from_sid_record(record):
        return [sid['compound'][-1]['id']['id']['cid'] for sid in record]
    
    def __cids_from_sid_records(self):
        self.cids = [Protein.cids_from_sid_record(df_sid_record) for df_sid_record in self.sid_records]
        
    @staticmethod
    def get_smiles_from_cids(cids):
        cids_str = ','.join([str(c) for c in cids])
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids_str}/property/CanonicalSMILES/JSON'
        response = requests.get(url)
        return response.json()['PropertyTable']['Properties'][0]['CanonicalSMILES']
    
    def __set_smiles(self):
        for cids in self.cids:
            self.smiles.append(Protein.get_smiles_from_cids(cids))
            time.sleep(0.5)
            
    def __add_smiles_to_dfs(self):
        for df, smiles in zip(self.concise_assay_dfs, self.smiles):
            df[self.__smiles_col_name] = smiles
    
    
if __name__ == '__main__':
    uniprot_id = 'P50129'
    protein = Protein(uniprot_id)
    print(protein.df)
