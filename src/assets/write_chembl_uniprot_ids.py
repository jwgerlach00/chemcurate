import requests
from typing import List


def get_chembl_uniprot_ids() -> List[str]:
    res = requests.get('https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_uniprot_mapping.txt')
    uniprot_ids = []
    for i, line in enumerate(res.text.splitlines()):
        if i == 0: # Header, skip
            chembl_version_num = line.split(' ')[1].split('_')[1]
        else:
            uniprot_ids.append(line.split('\t')[0])
    return uniprot_ids, chembl_version_num


if __name__ == '__main__':
    import joblib
    uniprot_ids = get_chembl_uniprot_ids()
    joblib.dump(uniprot_ids, 'chembl_uniprot_ids.joblib')
