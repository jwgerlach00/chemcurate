import requests
import pandas as pd


class Chembl:
    url_stem = 'https://www.ebi.ac.uk'
    
    def __init__(self, uniprot_ids:list):
        target_chembl_ids = [Chembl.uniprot_2_chembl_target_id(uniprot_id) for uniprot_id in uniprot_ids]
        self.df = Chembl.get_assays_from_target_chembl_ids(target_chembl_ids)
        structures = Chembl.get_structures_from_mol_chembl_ids(self.df['molecule_chembl_id'].tolist())
        self.df = structures.merge(self.df, on='molecule_chembl_id')
        
    def __str__(self):
        return str(self.df)
    
    def __len__(self) -> int:
        return len(self.df)
    
    @staticmethod
    def uniprot_2_chembl_target_id(uniprot_id:str):
        url = f'{Chembl.url_stem}/chembl/api/data/chembl_id_lookup/search.json?q={uniprot_id}'
        response = requests.get(url)
        return response.json()['chembl_id_lookups'][0]['chembl_id']
    
    @staticmethod
    def get_assays_from_target_chembl_ids(target_chembl_ids:list, limit:int=100) -> pd.DataFrame:
        url = '{0}/chembl/api/data/activity.json?target_chembl_id__in={1}&assay_type=B&limit={2}'\
            .format(Chembl.url_stem, ','.join(target_chembl_ids), limit)
        response = requests.get(url)
        
        activities = []
        activities += response.json()['activities']
        
        # Loop through remaining entries
        while response.json()['page_meta']['next']:
            response = requests.get('{0}{1}'.format(Chembl.url_stem, response.json()['page_meta']['next']))
            activities += response.json()['activities']
        
        columns_to_include = ['target_chembl_id', 'target_organism', 'target_pref_name', 'molecule_chembl_id',
                            'molecule_pref_name', 'pchembl_value', 'standard_type', 'standard_relation',
                            'standard_value', 'standard_units', 'assay_chembl_id','document_chembl_id', 'src_id']
        
        return pd.DataFrame(activities, columns=columns_to_include)

    @staticmethod
    def get_structures_from_mol_chembl_ids(mol_chembl_ids:list, limit:int=100) -> pd.DataFrame:
        url = '{0}/chembl/api/data/molecule.json?molecule_chembl_id__in={1}&limit={2}'\
            .format(Chembl.url_stem, ','.join(mol_chembl_ids), limit)
            
        response = requests.get(url)
        
        molecules = []
        molecules += response.json()['molecules']
        
        # Loop through remaining entries
        while response.json()['page_meta']['next']:
            response = requests.get('{0}{1}'.format(Chembl.url_stem, response.json()['page_meta']['next']))
            molecules += response.json()['molecules']
        
        canonical_smiles = [m['molecule_structures']['canonical_smiles'] for m in molecules]
        mol_chembl_ids = [m['molecule_chembl_id'] for m in molecules]
        
        return pd.DataFrame({'molecule_chembl_id': mol_chembl_ids, 'canonical_smiles': canonical_smiles})


if __name__ == '__main__':
    # target_id = 'CHEMBL2490'
    target_id = 'P50129'
    chembl = Chembl([target_id])
    print(chembl.df)
    print(len(chembl))
