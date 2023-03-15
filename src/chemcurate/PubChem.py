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


class PubChem:
    url_stem = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
    
    def __init__(self, uniprot_ids:list) -> None:
        self.__smiles_col_name = 'SMILES'
        
        self.uniprot_ids = uniprot_ids
        
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
        
    def __str__(self) -> str:
        return str(self.df)
    
    def __len__(self) -> int:
        return len(self.df)
    
    @staticmethod
    def get_assay_ids(uniprot_ids:list) -> list:
        url = '{0}/bioassay/target/ProteinName/{1}/aids/JSON'.format(PubChem.url_stem, ','.join(uniprot_ids))
        response = requests.get(url)
        return response.json()['IdentifierList']['AID']
    
    def __get_assay_ids(self) -> list:
        return PubChem.get_assay_ids(self.uniprot_ids)
    
    @staticmethod
    def get_concise_assay_dfs(assay_id:Union[int, str]) -> pd.DataFrame:
        url = f'{PubChem.url_stem}/bioassay/aid/{assay_id}/concise/JSON'
        response = requests.get(url)
        json = response.json()['Table']
        rows = [row['Cell'] for row in json['Row']]
        return pd.DataFrame(rows, columns=json['Columns']['Column'])

    def __set_concise_assay_dfs(self) -> None:
        for aid in self.assay_ids:
            self.concise_assay_dfs.append(PubChem.get_concise_assay_dfs(aid))
            time.sleep(0.5)
            
    @staticmethod
    def get_sid_records(sids) -> dict:
        sids_str = ','.join(sids)
        url = f'{PubChem.url_stem}/substance/sid/{sids_str}/record/JSON'
        response = requests.get(url)
        return response.json()['PC_Substances']
    
    def __set_sid_records(self) -> None:
        for df in self.concise_assay_dfs:
            self.sid_records.append(PubChem.get_sid_records(df['SID'].tolist()))
            time.sleep(0.5)
    
    @staticmethod
    def cids_from_sid_record(record) -> list:
        return [sid['compound'][-1]['id']['id']['cid'] for sid in record]
    
    def __cids_from_sid_records(self) -> None:
        self.cids = [PubChem.cids_from_sid_record(df_sid_record) for df_sid_record in self.sid_records]
        
    @staticmethod
    def get_smiles_from_cids(cids) -> list:
        cids_str = ','.join([str(c) for c in cids])
        url = f'{PubChem.url_stem}/compound/cid/{cids_str}/property/CanonicalSMILES/JSON'
        response = requests.get(url)
        return response.json()['PropertyTable']['Properties'][0]['CanonicalSMILES']
    
    def __set_smiles(self) -> None:
        for cids in self.cids:
            self.smiles.append(PubChem.get_smiles_from_cids(cids))
            time.sleep(0.5)
            
    def __add_smiles_to_dfs(self) -> None:
        for df, smiles in zip(self.concise_assay_dfs, self.smiles):
            df[self.__smiles_col_name] = smiles
    
    
if __name__ == '__main__':
    # uniprot_ids = ['P22303'] #['P50129', 'P10100']
    uniprot_ids = ['P50129']
    pc = PubChem(uniprot_ids)
    print(pc)
    print(len(pc))
