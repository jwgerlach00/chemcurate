import psycopg2
import pandas as pd


class PubChemQuery:
    def __init__(self):
        ########## Connect to DB ##########
        host = 'localhost'
        database = 'pubchem'
        user = 'postgres'
        password = 'Coll@bor@tions2020'
        # port = '5432'
        self.connection = psycopg2.connect(host=host, database=database, user=user, password=password)
        self.cursor = self.connection.cursor()
        
    def get_bioassays_from_uniprot_id(self, uniprot_id:str):
        self.cursor.execute(
            # NOTE: Assay data in index 0, bioassay_id is index 1, protein_accession is index 2, uniprot_id is index 3
            f'''
                SELECT json_array_elements(assay_data), bioassay_id, bioassay.protein_accession, test.uniprot_id
                FROM bioassay
                JOIN (
                    SELECT *
                    FROM target
                    WHERE target.uniprot_id = '{uniprot_id}'
                ) AS test ON bioassay.protein_accession = test.protein_accession;
            '''
        )
        results = self.cursor.fetchall()
        
        # TODO: user interface
        # TODO: SEPERATE BIOASSAY DATAFRAMES
        
        # First element: assay data (JSON) --> convert to dataframe
        results_df = pd.DataFrame.from_records(
            [x[0] for x in results] # index w/ 0 bc there is a list of tuples for some reason
        )
        
        # Add other fields to dataframe
        results_df['bioassay_id'] = [x[1] for x in results]
        results_df['protein_accession'] = [x[2] for x in results]
        results_df['uniprot_id'] = [x[3] for x in results]
        
        # Get SMILES from substance IDs
        substance_ids = tuple(results_df['sid'])
        self.cursor.execute(
            '''
                SELECT *
                FROM substance
                WHERE substance_id IN %s
            ''',
            (substance_ids,)
        )
        
        # Merges SMILES on substance IDs
        results_df = pd.merge(pd.DataFrame(self.cursor.fetchall(), columns=['sid', 'SMILES']), results_df, on='sid')
        return results_df
        

if __name__ == '__main__':
    pubchem_query = PubChemQuery()
    # pubchem_query.get_bioassays_from_uniprot_id('Q92952')
    pubchem_query.get_bioassays_from_uniprot_id('Q6W5P4')
