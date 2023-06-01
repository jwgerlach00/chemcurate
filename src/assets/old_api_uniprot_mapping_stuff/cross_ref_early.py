import yaml
import joblib
# from multiprocessing import Pool, cpu_count
import os

def process_organism(organism):
    print(organism)

    common_name = full_uniprot_mapping[organism]['common_name']
    uniprot_ids = [id for li in full_uniprot_mapping[organism]['protein'].values() for id in li]

    protein_dict = {}
    for uniprot_id in uniprot_ids:
        booly = (uniprot_id in chembl_uniprot_ids, uniprot_id in pubchem_uniprot_ids)
        if sum(booly):
            protein_dict[uniprot_id] = {
                'in_chembl': booly[0],
                'in_pubchem': booly[1]
            }

    if protein_dict:
        with open('cross_referenced_uniprot_mapping.yaml', 'a') as f:
            yaml.dump({
                organism: {
                    'common_name': common_name,
                    'protein': protein_dict
                }
            }, f, Dumper=yaml.SafeDumper)

if __name__ == '__main__':
    with open('uniprot_mapping.yaml') as f:
        full_uniprot_mapping = yaml.load(f, Loader=yaml.SafeLoader)

    chembl_uniprot_ids = set(joblib.load('chembl_uniprot_ids.joblib')[0]) # second element is chembl version
    pubchem_uniprot_ids = set(joblib.load('pubchem_uniprot_ids.joblib'))

    if os.path.exists('cross_referenced_uniprot_mapping.yaml'):
        os.remove('cross_referenced_uniprot_mapping.yaml')  

    [process_organism(organism) for organism in full_uniprot_mapping.keys()]
    # with Pool(cpu_count()) as p:
    #     p.map(process_organism, full_uniprot_mapping.keys())