from autochem.__asset_loader import uniprot_mapping
from autochem import PubChemQuery


if __name__ == '__main__':
    pubchem_query = PubChemQuery()
    uniprot_ids = pubchem_query.get_uniprot_ids_with_data()
    
    