import yaml
import joblib
# from multiprocessing import Pool, cpu_count
import os

def process_organism(self, organism:str):
    print(organism)
    common_name = self.full_uniprot_mapping[organism]['common_name']
    protein_names = [id for li in self.full_uniprot_mapping[organism]['protein'].keys() for id in li]
    uniprot_ids = [id for li in self.full_uniprot_mapping[organism]['protein'].values() for id in li]
    
    protein_dict = {}
    for protein_name in protein_names:
        uniprot_list = []
        for uniprot_id in uniprot_ids:
            booly = (uniprot_id in self.chembl_uniprot_ids, uniprot_id in self.pubchem_uniprot_ids)
            if sum(booly):
                uniprot_list.append({uniprot_id: {
                    'in_chembl': booly[0],
                    'in_pubchem': booly[1]
                }})
        if uniprot_list:
            protein_dict[protein_name]['common_name'] = common_name
            protein_dict[protein_name]['protein'] = uniprot_list

    if protein_dict:
        with open('cross_referenced_uniprot_mapping.yaml', 'a') as f:
            yaml.dump({organism: protein_dict}, f, Dumper=yaml.SafeDumper)

if __name__ == '__main__':
    with open('uniprot_mapping.yaml') as f:
        full_uniprot_mapping = yaml.load(f, Loader=yaml.SafeLoader)
       
    chembl_uniprot_ids = set(joblib.load('chembl_uniprot_ids.joblib')[0]) # second element is chembl version
    pubchem_uniprot_ids = set(joblib.load('pubchem_uniprot_ids.joblib'))
    
    if os.path.exists('cross_referenced_uniprot_mapping.yaml'):
        os.remove('cross_referenced_uniprot_mapping.yaml')  
    
    [process_organism(organism) for organism in full_uniprot_mapping.keys()]
