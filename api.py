import requests


def get_assay_ids(gene_symbol:str):
    # url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/caffeine/cids/TXT"
    # url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/490,1000/targets/ProteinGI,ProteinName,GeneID,GeneSymbol/XML'
    # url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/3193/record/CSV'
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/genesymbol/{gene_symbol}/aids/JSON'
    response = requests.get(url)
    return response.json()['IdentifierList']['AID']

def get_assay_records(assay_id:int):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{assay_id}/record/CSV'
    response = requests.get(url)
    return response.text


if __name__ == '__main__':
    ids = get_assay_ids('5-HT2A')
    out = get_assay_records(ids[2])
    print(out)