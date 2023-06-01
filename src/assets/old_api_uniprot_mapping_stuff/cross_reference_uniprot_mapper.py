import yaml
import joblib
# from multiprocessing import Pool, cpu_count
import os

def process_organism(organism:str):
    print(organism)
    protein_names = full_uniprot_mapping[organism]['protein'].keys()
    
    protein_dict = {}
    for protein_name in protein_names:
        uniprot_dict = {}
        for uniprot_id in full_uniprot_mapping[organism]['protein'][protein_name]:
            booly = (uniprot_id in chembl_uniprot_ids, uniprot_id in pubchem_uniprot_ids)
            if sum(booly):
                uniprot_dict[uniprot_id] =  {
                    'in_chembl': booly[0],
                    'in_pubchem': booly[1]
                }
        if uniprot_dict:
            protein_dict[protein_name] = uniprot_dict
            
    if protein_dict:
        return {
            'common_name': full_uniprot_mapping[organism]['common_name'],
            'protein': protein_dict
        }
    else:
        return None

    # if protein_dict:
    #     with open('cross_referenced_uniprot_mapping.yaml', 'a') as f:
    #         yaml.dump({organism: protein_dict}, f, Dumper=yaml.SafeDumper)

if __name__ == '__main__':
    with open('uniprot_mapping.yaml') as f:
        full_uniprot_mapping = yaml.load(f, Loader=yaml.SafeLoader)
       
    chembl_uniprot_ids = set(joblib.load('chembl_uniprot_ids.joblib')[0]) # second element is chembl version
    pubchem_uniprot_ids = set(joblib.load('pubchem_uniprot_ids.joblib'))
    
    if os.path.exists('cross_referenced_uniprot_mapping.yaml'):
        os.remove('cross_referenced_uniprot_mapping.yaml')  
    
    for organism in full_uniprot_mapping.keys():
        x = process_organism(organism)
        if x:
            with open('cross_referenced_uniprot_mapping.yaml', 'a') as f:
                yaml.dump({organism: x}, f, Dumper=yaml.SafeDumper)
        else:
            print('ape')
    # [process_organism(organism) for organism in full_uniprot_mapping.keys()]
