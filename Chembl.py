import requests
import pandas as pd


class Chembl:
    def __init__(self, target_chembl_ids:list):
        self.df = Chembl.get_assays_from_target_chembl_ids(target_chembl_ids)
        structures = Chembl.get_structures_from_mol_chembl_ids(self.df['molecule_chembl_id'].tolist())
        self.df = structures.merge(self.df, on='molecule_chembl_id')
        
    def __str__(self):
        return str(self.df)
    
    def __len__(self) -> int:
        return len(self.df)
    
    @staticmethod
    def get_assays_from_target_chembl_ids(target_chembl_ids:list, limit:int=100) -> pd.DataFrame:
        url = 'https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id__in={0}&assay_type=B&limit={1}'\
            .format(','.join(target_chembl_ids), limit)
        response = requests.get(url)
        
        activities = []
        activities += response.json()['activities']
        
        # Loop through remaining entries
        while response.json()['page_meta']['next']:
            response = requests.get('https://www.ebi.ac.uk{0}'.format(response.json()['page_meta']['next']))
            activities += response.json()['activities']
        
        columns_to_include = ['target_chembl_id', 'target_organism', 'target_pref_name', 'molecule_chembl_id',
                            'molecule_pref_name', 'pchembl_value', 'standard_type', 'standard_relation',
                            'standard_value', 'standard_units', 'assay_chembl_id','document_chembl_id', 'src_id']
        
        return pd.DataFrame(activities, columns=columns_to_include)

    @staticmethod
    def get_structures_from_mol_chembl_ids(mol_chembl_ids:list, limit:int=100) -> pd.DataFrame:
        url = 'https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_chembl_id__in={0}&limit={1}'\
            .format(','.join(mol_chembl_ids), limit)
            
        response = requests.get(url)
        
        molecules = []
        molecules += response.json()['molecules']
        
        # Loop through remaining entries
        while response.json()['page_meta']['next']:
            response = requests.get('https://www.ebi.ac.uk{0}'.format(response.json()['page_meta']['next']))
            molecules += response.json()['molecules']
        
        canonical_smiles = [m['molecule_structures']['canonical_smiles'] for m in molecules]
        mol_chembl_ids = [m['molecule_chembl_id'] for m in molecules]
        
        return pd.DataFrame({'molecule_chembl_id': mol_chembl_ids, 'canonical_smiles': canonical_smiles})

if __name__ == '__main__':
    target_id = 'CHEMBL2490'
    chembl = Chembl([target_id])
    print(chembl.df)
    print(len(chembl))