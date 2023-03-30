import requests
import pandas as pd
from typing import Union, List
import time

from chemcurate import __Base


# def get_full_bioassay_record(assay_id:Union[int, str]):
#     url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/bioassay/aid/{assay_id}/record/CSV'
#     response = requests.get(url)
#     return response.text

# def get_assay_summary(assay_id:Union[int, str]):
#     url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/bioassay/aid/{assay_id}/summary/JSON'
#     response = requests.get(url)
#     return response.text


class PubChem(__Base):
    _url_stem = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
    
    def __init__(self, uniprot_ids:list) -> None:
        super(PubChem, self).__init__()
        
        self.__smiles_col_name = 'SMILES'
        
        self.cids = []
        self.smiles = []
        
        assay_ids = PubChem.get_assay_ids(uniprot_ids)
        concise_assay_dfs = PubChem.get_concise_assay_df(assay_ids)
        sid_records = [PubChem.get_sid_records(df['SID'].tolist()) for df in concise_assay_dfs]
        cids = [PubChem.cids_from_sid_record(df_sid_record) for df_sid_record in sid_records]
    
        for cids in self.cids:
            self.smiles.append(PubChem.get_smiles_from_cids(cids))
            time.sleep(0.5)
        
        
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
    def __sleep() -> None:
        time.sleep(0.5)
    
    @staticmethod
    def get_assay_ids(uniprot_ids:List[str]) -> List[str]:
        aids = []
        for batch in PubChem.batch(uniprot_ids, PubChem._batch_size):
            url = '{0}/bioassay/target/ProteinName/{1}/aids/JSON'.format(PubChem._url_stem, ','.join(batch))
            aids += requests.get(url).json()['IdentifierList']['AID']
            PubChem.__sleep()
        return aids
    
    @staticmethod
    def get_concise_assay_dfs(assay_ids:List[Union[int, str]]) -> pd.DataFrame:
        dfs = []
        for batch in PubChem.batch(assay_ids, PubChem._batch_size):
            url = '{0}/bioassay/aid/{1}/concise/JSON'.format(PubChem._url_stem, ','.join(batch))
            response = requests.get(url)
            json = response.json()['Table']
            rows = [row['Cell'] for row in json['Row']]
            dfs.append(pd.DataFrame(rows, columns=json['Columns']['Column']))
            PubChem.__sleep()
        return dfs
            
    @staticmethod
    def get_sid_records(sids:List[Union[int, str]]) -> List[dict]:        
        records = []
        for batch in PubChem.batch(sids, PubChem._batch_size):
            url = '{0}/substance/sid/{1}/record/JSON'.format(PubChem._url_stem, ','.join(batch))
            records += requests.get(url).json()['PC_Substances']
            PubChem.__sleep()
        return records
    
    @staticmethod
    def cids_from_sid_record(sid_record:dict) -> list:
        return [sid['compound'][-1]['id']['id']['cid'] for sid in sid_record]
        
    @staticmethod
    def get_smiles_from_cid(cids:List[str]) -> List[str]:
        smiles = []
        for batch in PubChem.batch(cids, PubChem._batch_size):
            url = '{0}/compound/cid/{cid}/property/CanonicalSMILES/JSON'.format(PubChem._url_stem, ','.join(batch))
            smiles += requests.get(url).json()['PropertyTable']['Properties'][0]['CanonicalSMILES']
            PubChem.__sleep()
        return smiles
    
    def __set_smiles(self) -> None:
        for cids in self.cids:
            self.smiles.append(PubChem.get_smiles_from_cids(cids))
            time.sleep(0.5)
            
    def __add_smiles_to_dfs(self) -> None:
        for df, smiles in zip(self.concise_assay_dfs, self.smiles):
            df[self.__smiles_col_name] = smiles
    
    
if __name__ == '__main__':
    # uniprot_ids = ['P22303'] #['P50129', 'P10100']
    uniprot_ids = ['P50129']#, 'Q13936']
    pc = PubChem(uniprot_ids)
    print(pc)
    print(len(pc))
