import yaml
import autochem
import os
import tqdm
import itertools

if __name__ == '__main__':
    yaml.Dumper.ignore_aliases = lambda *args : True
    
    with open('uniprot_mapping.yaml') as f:
        uniprot_mapping = yaml.load(f, Loader=yaml.SafeLoader)
    
    with open('chembl_uniprot_ids.yaml') as f:
        chembl_uniprot_ids = yaml.load(f, Loader=yaml.SafeLoader)['uniprot_ids']

    if os.path.exists('cross_uniprot_chembl.yaml'):
        os.remove('cross_uniprot_chembl.yaml')
    
    for organism in uniprot_mapping.keys():
        print(organism)
        
        common_name = uniprot_mapping[organism]['common_name']
        
        protein_names, uniprot_ids = uniprot_mapping[organism]['protein'].keys(), \
            itertools.chain.from_iterable(uniprot_mapping[organism]['protein'].values())
        
        out_organism_dict = {}
        out_uniprot_ids = []
        for protein_name in protein_names:
            for uniprot_id in tqdm.tqdm(uniprot_ids):
                in_chembl = uniprot_id in chembl_uniprot_ids
                if in_chembl:
                    # organism_out_dict[uniprot_id] = {'in_chembl': in_chembl, 'in_pubchem': None}
                    out_uniprot_ids.append(uniprot_id)
            if out_uniprot_ids:
                out_organism_dict[protein_name] = out_uniprot_ids
            
        # Dump once for every protein
        if out_organism_dict:
            with open('cross_uniprot_chembl.yaml', 'a') as f:
                yaml.dump({
                    organism: {
                        'common_name': common_name,
                        'protein': out_organism_dict
                    }
                }, f)
