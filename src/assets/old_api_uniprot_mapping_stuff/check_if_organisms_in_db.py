import yaml
from yaml import SafeLoader
import autochem
import os
import tqdm


if __name__ == '__main__':
    with open('uniprot_mapping_filtered.yaml') as f:
        uniprot_mapping = yaml.load(f, Loader=SafeLoader)

    if os.path.exists('organisms_in_db.yaml'):
        os.remove('organisms_in_db.yaml')
    
    for organism in uniprot_mapping.keys():
        print(organism)
        
        uniprot_ids = []
        for li in uniprot_mapping['Homo sapiens']['protein'].values():
            uniprot_ids += li

        organism_out_dict = {}
        for uniprot_id in tqdm.tqdm(uniprot_ids):
            try:
                in_chembl = bool(autochem.Chembl.uniprot_2_chembl_target_id(uniprot_id))
                in_pubchem = bool(autochem.PubChem.get_assay_ids([uniprot_id]))
            except:
                in_chembl = False
                in_pubchem = False
            
            organism_out_dict[uniprot_id] = {'in_chembl': in_chembl, 'in_pubchem': in_pubchem}
            
        # Dump once for every organism
        with open('organisms_in_db.yaml', 'a') as f:
            yaml.dump({organism: organism_out_dict}, f)

    