import requests
from requests.adapters import HTTPAdapter, Retry
import re
from copy import deepcopy
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pprint import pprint


class UniprotMapper:
    # Sourced from https://www.uniprot.org/help/api_queries
    _re_next_link = re.compile(r'<(.+)>; rel="next"')
    _session = requests.Session()
    _session.mount("https://", HTTPAdapter(
        max_retries=Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    ))

    _full_uniprot_url = 'https://rest.uniprot.org/uniprotkb/search?fields=accession,protein_name,organism_name&format=json&query=(reviewed:true)&size=500'

    def __init__(self) -> None:
        pass

    @staticmethod
    def get_full_uniprot_map(verbose:bool=True) -> dict:
        uniprot_mapping = {}
        for batch, total in UniprotMapper.get_batch(UniprotMapper._full_uniprot_url):
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
                            
                protein_common_name = result['proteinDescription']['recommendedName']['fullName']['value']
                
                uniprot_id = result['primaryAccession']
                uniprot_mapping[organism_sci_name]['protein'][uniprot_id] = {'common_name': protein_common_name}
            
            if verbose:
                pprint(uniprot_mapping)
                current_len = sum([len(uniprot_mapping[key]['protein']) for key in uniprot_mapping.keys()])
                print(f'{current_len} / {total}')
        
        return uniprot_mapping
    
    @staticmethod
    def get_batch(batch_url:str):
        '''Sourced from: https://www.uniprot.org/help/api_queries'''
        while batch_url:
            response = UniprotMapper._session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = UniprotMapper.get_next_link(response.headers)
            
    @staticmethod
    def get_next_link(headers):
        '''Sourced from: https://www.uniprot.org/help/api_queries'''
        if 'Link' in headers:
            match = UniprotMapper._re_next_link.match(headers["Link"])
            if match:
                return match.group(1)


if __name__ == '__main__':
    uniprot_mapper = UniprotMapper()
    uniprot_mapper.get_full_uniprot_map()
