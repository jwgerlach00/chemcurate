import psycopg2
import requests
from requests.adapters import HTTPAdapter, Retry
import re
from tqdm import tqdm
from .__ABCChemDB import __ABCChemDB
from math import ceil


class KnownLengthIterator:
    def __init__(self, it, length):
        self.it = it
        self.length = int(length)

    def __len__(self):
        return self.length

    def __iter__(self):
        yield from self.it


class UniProtDB(__ABCChemDB):
    def __init__(self) -> None:
        ########## Connect to DB ##########
        host = 'localhost'
        database = 'uniprot'
        user = 'postgres'
        password = 'Coll@bor@tions2020'
        # port = '5432'
        self.connection = psycopg2.connect(host=host, database=database, user=user, password=password)
        self.cursor = self.connection.cursor()
        
        # Sourced from https://www.uniprot.org/help/api_queries
        self.__re_next_link = re.compile(r'<(.+)>; rel="next"')
        self.__session = requests.Session()
        
        retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
        self.__session.mount("https://", HTTPAdapter(max_retries=retries))
        
        self.__url = 'https://rest.uniprot.org/uniprotkb/search?fields=accession,protein_name,organism_name&format=json&query=(reviewed:true)&size=500'
        
    def build(self) -> None:
        # Clear table
        self.cursor.execute('DELETE FROM organism_uniprot')
        self.connection.commit()
        
        n_batches = ceil(int(next(self._get_batch(self.__url))[1])/500)
        # generator = KnownLengthIterator(self._get_batch(self.__url), n_batches)
        
        for batch, _ in tqdm(self._get_batch(self.__url), total=n_batches):
            results = batch.json()['results']
            
            for result in results:
                uniprot_id = result['primaryAccession']
                protein_name = result['proteinDescription']['recommendedName']['fullName']['value']
                organism_sci_name = result['organism']['scientificName']
                organism_common_name = result['organism']['commonName'] if 'commonName' in result['organism'] else None
                
                self.cursor.execute('INSERT INTO organism_uniprot (uniprot_id, protein_name, organism_sci_name, \
                    organism_common_name) VALUES (%s, %s, %s, %s)', (uniprot_id, protein_name, organism_sci_name,
                                                                     organism_common_name))
                self.connection.commit()

    def _get_batch(self, batch_url):
        '''Sourced from: https://www.uniprot.org/help/api_queries'''
        while batch_url:
            response = self.__session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = self._get_next_link(response.headers)
            
    def _get_next_link(self, headers):
        '''Sourced from: https://www.uniprot.org/help/api_queries'''
        if "Link" in headers:
            match = self.__re_next_link.match(headers["Link"])
            if match:
                return match.group(1)
