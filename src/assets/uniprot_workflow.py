from chemcurate import UniprotMapper
import joblib
import yaml


if __name__ == '__main__':
    chembl_uniprot_ids = joblib.load('chembl_uniprot_ids.joblib')[0]
    pubchem_uniprot_ids = joblib.load('pubchem_uniprot_ids.joblib')
    out_dict = UniprotMapper.get_processed_uniprot_map(chembl_uniprot_ids, pubchem_uniprot_ids, verbose=True, save_path='ape.yaml')
    with open('test.yaml', 'w') as f:
        yaml.dump(out_dict, f, Dumper=yaml.SafeDumper)
    
    