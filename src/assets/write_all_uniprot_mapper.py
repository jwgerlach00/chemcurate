# import requests#, Retry

# # retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
# # session = requests.Session()
# # session.mount("https://", HTTPAdapter(max_retries=retries))


# # def get_next_link(headers):
# #     if "Link" in headers:
# #         match = re_next_link.match(headers["Link"])
# #         if match:
# #             return match.group(1)

# # def get_batch(batch_url):
# #     while batch_url:
# #         response = session.get(batch_url)
# #         response.raise_for_status()
# #         total = response.headers["x-total-results"]
# #         yield response, total
# #         batch_url = get_next_link(response.headers)


# url = 'https://rest.uniprot.org/uniprotkb/search?query=reviewed:true&size=500&format=json'

# response = requests.get(url)
# response_json = response.json()
# response_headers = response.headers




# from pprint import pprint
# pprint([result['primaryAccession'] for result in response_json['results']])
# pprint(len(response_json['results']))
# pprint(response_headers)

import requests
from requests.adapters import HTTPAdapter, Retry
import re
from pprint import pprint
import yaml


# Sourced from https://www.uniprot.org/help/api_queries
re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def get_next_link(headers):
    '''Sourced from: https://www.uniprot.org/help/api_queries'''
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    '''Sourced from: https://www.uniprot.org/help/api_queries'''
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)


if __name__ == '__main__':
    url = 'https://rest.uniprot.org/uniprotkb/search?fields=accession,protein_name,organism_name&format=json&query= \
        (reviewed:true)&size=500'

    uniprot_mapping = {}
    for batch, total in get_batch(url):
        results = batch.json()['results']
        
        for result in results:
            organism_sci_name = result['organism']['scientificName']
            if organism_sci_name not in uniprot_mapping:
                uniprot_mapping[organism_sci_name] = {}
                uniprot_mapping[organism_sci_name]['protein'] = {}    
                
            
            if 'commonName' in result['organism']:
                organism_common_name = result['organism']['commonName']
                uniprot_mapping[organism_sci_name]['common_name'] = organism_common_name
            else:
                uniprot_mapping[organism_sci_name]['common_name'] = None
                        
            protein_name = result['proteinDescription']['recommendedName']['fullName']['value']
            # print(uniprot_mapping)
            if protein_name not in uniprot_mapping[organism_sci_name]['protein']:
                uniprot_mapping[organism_sci_name]['protein'][protein_name] = [result['primaryAccession']]  # Turn into list to allow for multiple mappings
            else:
                uniprot_mapping[organism_sci_name]['protein'][protein_name].append(result['primaryAccession'])
        
        current_len = sum([len(uniprot_mapping[key]['protein']) for key in uniprot_mapping.keys()])
        print(f'{current_len} / {total}')
        # break
    
    with open('uniprot_mapping.yaml', 'w') as f:
        yaml.dump(uniprot_mapping, f)
