import requests
from autochem import PubChem
import joblib


def get_pubchem_uniprot_ids():
    url = lambda page: f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/annotations/heading/JSON/?source=UniProt& \
        heading_type=Protein&heading=UniProt%20ID&page={page}'
    
    out = []
    page = 1
    total_pages = 2 # Some arbitrary number greater than page to start
    while page < total_pages:
        print(page)
        res_json = requests.get(url(page)).json()
        if page == 1:
            total_pages = res_json['Annotations']['TotalPages']
        for i in res_json['Annotations']['Annotation']:
            out.extend(i['LinkedRecords']['ProteinAccession'])
        page += 1
        
        PubChem._sleep()

    return out
    
if __name__ == '__main__':
    uniprot_ids = get_pubchem_uniprot_ids()
    print(uniprot_ids)
    joblib.dump(uniprot_ids, 'pubchem_uniprot_ids.joblib')
