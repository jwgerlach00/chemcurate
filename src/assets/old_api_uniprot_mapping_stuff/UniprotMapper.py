import requests
from requests.adapters import HTTPAdapter, Retry
import re
from chemcurate.api_tools import Chembl, PubChem
from copy import deepcopy
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed


class UniprotMapper:
    # Sourced from https://www.uniprot.org/help/api_queries
    _re_next_link = re.compile(r'<(.+)>; rel="next"')
    _session = requests.Session()
    _session.mount("https://", HTTPAdapter(
        max_retries=Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    ))

    _full_uniprot_url = 'https://rest.uniprot.org/uniprotkb/search?fields=accession,protein_name,organism_name&format=json&query=(reviewed:true)&size=500'

    def __init__(self, full_uniprot_map=None, chembl_uniprot_ids=None, pubchem_uniprot_ids=None):
        self.full_uniprot_map = deepcopy(full_uniprot_map) if full_uniprot_map else UniprotMapper.get_full_uniprot_map()
        self.chembl_uniprot_ids = deepcopy(chembl_uniprot_ids) if chembl_uniprot_ids else \
            set(Chembl.get_db_uniprot_ids()[0]) # second element is chembl version
        self.pubchem_uniprot_ids = deepcopy(pubchem_uniprot_ids) if pubchem_uniprot_ids else \
            set(PubChem.get_db_uniprot_ids())
        
        self.map = {organism: self.process_organism(organism) for organism in self.full_uniprot_map.keys()}
        del self.full_uniprot_map
        del self.chembl_uniprot_ids
        del self.pubchem_uniprot_ids
        
        return self.map
        
    @staticmethod
    def get_next_link(headers):
        '''Sourced from: https://www.uniprot.org/help/api_queries'''
        if 'Link' in headers:
            match = UniprotMapper._re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

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
    def get_full_uniprot_map(verbose:bool=False):
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
                print(uniprot_mapping)
                current_len = sum([len(uniprot_mapping[key]['protein']) for key in uniprot_mapping.keys()])
                print(f'{current_len} / {total}')
        
        return uniprot_mapping

    @staticmethod
    def process_uniprot_batch(batch):
        results = batch.json()['results']
        organism_sci_names = [result['organism']['scientificName'] for result in results]
        organism_common_names = [result['organism'].get('commonName') for result in results]
        protein_common_names = [result['proteinDescription']['recommendedName']['fullName']['value'] for result in results]
        uniprot_ids = [result['primaryAccession'] for result in results]
        batch_results = []
        for organism_sci_name, organism_common_name, protein_common_name, uniprot_id in zip(organism_sci_names, organism_common_names, protein_common_names, uniprot_ids):
            batch_results.append((uniprot_id, organism_sci_name, organism_common_name, protein_common_name))
        return batch_results

    @staticmethod
    def get_processed_uniprot_map(chembl_uniprot_ids, pubchem_uniprot_ids, verbose=False, save_path=None, num_threads=4):
        uniprot_mapping = {
            'chembl': {},
            'pubchem': {},
        }
        
        chembl_uniprot_ids_set = set(chembl_uniprot_ids)
        pubchem_uniprot_ids_set = set(pubchem_uniprot_ids)
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            tasks = []
            for batch, total in UniprotMapper.get_batch(UniprotMapper._full_uniprot_url):
                task = executor.submit(UniprotMapper.process_uniprot_batch, batch)
                tasks.append(task)
                print('run1')
            for task in as_completed(tasks):
                print('run2')
                batch_results = task.result()
                for uniprot_id, organism_sci_name, organism_common_name, protein_common_name in batch_results:
                    if uniprot_id in chembl_uniprot_ids_set:
                        db_str = 'chembl'
                    elif uniprot_id in pubchem_uniprot_ids_set:
                        db_str = 'pubchem'
                    else:
                        continue
                    if organism_sci_name not in uniprot_mapping[db_str]:
                        if verbose:
                            print(organism_sci_name)
                        uniprot_mapping[db_str][organism_sci_name] = {'common_name': organism_common_name, 'protein': {}}
                    protein_dict = uniprot_mapping[db_str][organism_sci_name]['protein']
                    if uniprot_id not in protein_dict:
                        protein_dict[uniprot_id] = {'common_name': protein_common_name}
                # if verbose:
                #     print(uniprot_mpaaing)
                #     current_len = sum([len(uniprot_mapping[key]['protein']) for key in uniprot_mapping.keys()])
                #     print(f'{current_len} / {total}')
        return uniprot_mapping


    
    # @staticmethod
    # def get_processed_uniprot_map(chembl_uniprot_ids, pubchem_uniprot_ids, verbose:bool=False, save_path:str=None):
    #     uniprot_mapping = {
    #         'chembl': {},
    #         'pubchem': {},
    #     }
    #     for batch,  total in UniprotMapper.get_batch(UniprotMapper._full_uniprot_url):
    #         results = batch.json()['results']
            
    #         for result in results:
    #             organism_sci_name = result['organism']['scientificName']
    #             # if organism_sci_name not in uniprot_mapping:
    #             #     uniprot_mapping[organism_sci_name] = {}
    #             #     uniprot_mapping[organism_sci_name]['protein'] = {}    
                    
    #             # if 'commonName' in result['organism']:
    #             organism_common_name = result['organism']['commonName'] if 'commonName' in result['organism'] else None
    #                 # uniprot_mapping[organism_sci_name]['common_name'] = organism_common_name
    #             # else:
    #             #     uniprot_mapping[organism_sci_name]['common_name'] = None
                            
    #             protein_common_name = result['proteinDescription']['recommendedName']['fullName']['value']
                
    #             uniprot_id = result['primaryAccession']
                
    #             def check_organism_in(organism_sci_name, organism_common_name, db_str):
    #                 if organism_sci_name not in uniprot_mapping[db_str]:
    #                     print(organism_sci_name)
    #                     uniprot_mapping[db_str][organism_sci_name] = {}
    #                     uniprot_mapping[db_str][organism_sci_name]['common_name'] = organism_common_name
    #                     uniprot_mapping[db_str][organism_sci_name]['protein'] = {}
                
    #             if uniprot_id in chembl_uniprot_ids:
    #                 check_organism_in(organism_sci_name, organism_common_name, 'chembl')
    #                 uniprot_mapping['chembl'][organism_sci_name]['protein'][uniprot_id] = {'common_name': protein_common_name}
    #             if uniprot_id in pubchem_uniprot_ids:
    #                 check_organism_in(organism_sci_name, organism_common_name, 'pubchem')
    #                 uniprot_mapping['pubchem'][organism_sci_name]['protein'][uniprot_id] = {'common_name': protein_common_name}
                
    #             # uniprot_mapping[organism_sci_name]['protein'][uniprot_id] = {'common_name': protein_common_name}
    #             # print(uniprot_mapping.keys())
            
    #         # if verbose:
    #         #     print()
    #             # print(uniprot_mapping)
    #             # current_len = sum([len(uniprot_mapping[key]['protein']) for key in uniprot_mapping.keys()])
    #             # print(f'{current_len} / {total}')
                
    #     return uniprot_mapping

    def process_organism(self, organism:str, verbose:bool=False):
        if verbose:
            print(organism)

        protein_names = self.full_uniprot_mapping[organism]['protein'].keys()
        
        protein_dict = {}
        for protein_name in protein_names:
            uniprot_list = []
            for uniprot_id in self.full_uniprot_mapping[organism]['protein'][protein_name]:
                booly = (uniprot_id in self.chembl_uniprot_ids, uniprot_id in self.pubchem_uniprot_ids)
                if sum(booly):
                    uniprot_list.append({uniprot_id: {
                        'in_chembl': booly[0],
                        'in_pubchem': booly[1]
                    }})
            if uniprot_list:
                protein_dict[protein_name] = uniprot_list
                
        if protein_dict:
            return {
                'common_name': self.full_uniprot_mapping[organism]['common_name'],
                'protein': protein_dict
            }
        else:
            return None
