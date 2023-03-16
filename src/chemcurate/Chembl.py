import requests
import pandas as pd
from typing import Union, List, Generator, Iterable


class Chembl:
    _url_stem = 'https://www.ebi.ac.uk'
    _batch_size = 100 # NOTE: DONT CHANGE (Chembl API limit)
    
    def __init__(self, uniprot_ids:List[str]) -> None:
        target_chembl_ids = []
        for i in uniprot_ids:
            target_chembl_ids += Chembl.uniprot_2_chembl_target_id(i)

        target_chembl_ids = [i for i in target_chembl_ids if i] # Remove None values
        if target_chembl_ids:
            self.df = Chembl.get_assays_from_target_chembl_ids(target_chembl_ids)
            structures = Chembl.get_structures_from_mol_chembl_ids(self.df['molecule_chembl_id'].tolist())
            self.df = structures.merge(self.df, on='molecule_chembl_id')
        else:
            self.df = pd.DataFrame() # Empty dataframe if there is no data
        
    def __str__(self) -> str:
        return str(self.df)
    
    def __len__(self) -> int:
        return len(self.df)
    
    @staticmethod
    def batch(iterable:Iterable, n:int=1) -> Generator:
        l = len(iterable)
        for ndx in range(0, l, n):
            yield iterable[ndx:min(ndx + n, l)]
    
    @staticmethod
    def uniprot_2_chembl_target_id(uniprot_id:str) -> Union[List[str], None]:
        url = f'{Chembl._url_stem}/chembl/api/data/chembl_id_lookup/search.json?q={uniprot_id}'
        response_json = requests.get(url).json()
        if response_json['chembl_id_lookups']:
            return [target['chembl_id'] for target in response_json['chembl_id_lookups']]
        else:
            return None
        
    @staticmethod
    def get_assays_from_target_chembl_ids(target_chembl_ids:List[str]) -> pd.DataFrame:
        activities = []
        # Batch in case targets are greater than Chembl._batch_size
        for x in Chembl.batch(target_chembl_ids, Chembl._batch_size):
            url = '{0}/chembl/api/data/activity.json?target_chembl_id__in={1}&assay_type=B&limit={2}'\
                .format(Chembl._url_stem, ','.join(x), Chembl._batch_size)
            
            response_json = requests.get(url).json()
            activities += response_json['activities']
            
            # Loop through remaining entries
            while response_json['page_meta']['next']:
                response_json = requests.get('{0}{1}'.format(Chembl._url_stem,
                                                             response_json['page_meta']['next'])).json()
                activities += response_json['activities']
        
        columns_to_include = ['target_chembl_id', 'target_organism', 'target_pref_name', 'molecule_chembl_id',
                            'molecule_pref_name', 'pchembl_value', 'standard_type', 'standard_relation',
                            'standard_value', 'standard_units', 'assay_chembl_id','document_chembl_id', 'src_id']
        
        return pd.DataFrame(activities, columns=columns_to_include)
    
    @staticmethod
    def get_structures_from_mol_chembl_ids(mol_chembl_ids:List[str]) -> pd.DataFrame:
        unique_mol_chembl_ids = list(set(mol_chembl_ids))
        
        molecules = []
        for x in Chembl.batch(unique_mol_chembl_ids, Chembl._batch_size):
            url = '{0}/chembl/api/data/molecule.json?molecule_chembl_id__in={1}&limit={2}'\
                .format(Chembl._url_stem, ','.join(x), Chembl._batch_size)
            response_json = requests.get(url).json()
            molecules += response_json['molecules']
            # NOTE: No next page because retults = limit (one-to-one)
        
        mcid = []
        cs = []
        [(mcid.append(m['molecule_chembl_id']),
          cs.append(m['molecule_structures']['canonical_smiles']))
         for m in molecules]
            
        return pd.DataFrame({
            'molecule_chembl_id': mcid,
            'canonical_smiles': cs
        })


if __name__ == '__main__':
    target_id_1 = 'P50129' # chembl_id = CHEMBL2490
    target_id_2 = 'Q13936'
    chembl = Chembl([target_id_1, target_id_2])
    print(chembl.df)
    print(len(chembl))
