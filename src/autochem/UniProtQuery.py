import psycopg2
from typing import Tuple
import pandas as pd


class UniProtQuery:
    def __init__(self) -> None:
        ########## Connect to DB ##########
        host = 'localhost'
        database = 'uniprot'
        user = 'postgres'
        password = 'Coll@bor@tions2020'
        # port = '5432'
        self.connection = psycopg2.connect(host=host, database=database, user=user, password=password)
        self.cursor = self.connection.cursor()
        
    def get_all_data_for_uniprot_ids(self, uniprot_ids:Tuple[str]) -> pd.DataFrame:
        self.cursor.execute(
            '''
                SELECT *
                FROM organism_uniprot
                WHERE uniprot_id IN %s
            ''',
            (uniprot_ids,)
        )
        return pd.DataFrame(self.cursor.fetchall(), columns=['uniprot_id', 'protein_name', 'organism_sci_name',
                                                             'organism_common_name'])
