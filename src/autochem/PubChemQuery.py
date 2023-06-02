import psycopg2
import pandas as pd
from typing import List, Dict
from autochem.UniProtQuery import UniProtQuery


class PubChemQuery:
    def __init__(self) -> None:
        ########## Connect to DB ##########
        host = 'localhost'
        database = 'pubchem'
        user = 'postgres'
        password = 'Coll@bor@tions2020'
        # port = '5432'
        self.connection = psycopg2.connect(host=host, database=database, user=user, password=password)
        self.cursor = self.connection.cursor()
        
    def get_assay_data_from_bioassay_id(self, bioassay_id:str):
        self.cursor.execute(
            f'''
                SELECT json_array_elements(assay_data)
                FROM bioassay
                WHERE bioassay_id = '{bioassay_id}'
            '''
        )
        out = self.cursor.fetchall()
        from pprint import pprint
        pprint(out)
        
    def get_bioassays_from_protein_accession(self, protein_accession:str) -> pd.DataFrame:
        self.cursor.execute(
            # NOTE: Assay data in index 0, bioassay_id is index 1, protein_accession is index 2, uniprot_id is index 3
            f'''
                SELECT json_array_elements(assay_data), bioassay_id, protein_accession
                FROM bioassay
                WHERE protein_accession = '{protein_accession}' AND assay_data IS NOT NULL AND json_typeof(assay_data)='array'
            '''
        )
        results = self.cursor.fetchall()
        
        # First element: assay data (JSON) --> convert to dataframe
        results_df = pd.DataFrame.from_records(
            [x[0] for x in results] # index w/ 0 bc there is a list of tuples for some reason
        )
        
        # Add other fields to dataframe
        results_df['bioassay_id'] = [x[1] for x in results]
        results_df['protein_accession'] = [x[2] for x in results]
        
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
        
    def get_bioassays_from_uniprot_id(self, uniprot_id:str) -> pd.DataFrame:
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
    
    def get_represented_uniprot_ids(self) -> List[str]:
        self.cursor.execute(
            '''
                SELECT DISTINCT uniprot_id
                FROM target
            '''
        )
        return [x[0] for x in self.cursor.fetchall()]
    
    def get_uniprot_ids_with_data(self) -> List[str]:
        self.cursor.execute(
            '''
                SELECT DISTINCT target.uniprot_id
                FROM target
                JOIN bioassay
                ON target.protein_accession = bioassay.protein_accession
                WHERE bioassay.assay_data IS NOT NULL AND target.uniprot_id IS NOT NULL
            '''
        )
        return [x[0] for x in self.cursor.fetchall()]
    
    def get_protein_accessions_with_data(self) -> List[str]:
        self.cursor.execute(
            '''
                SELECT DISTINCT protein_accession
                FROM bioassay
                WHERE protein_accession IS NOT NULL AND protein_acces
            '''
        )
        return [x[0] for x in self.cursor.fetchall()]
    
    def get_bioassay_protein_accessions(self) -> List[str]:
        self.cursor.execute(
            '''
                SELECT DISTINCT protein_accession
                FROM bioassay
                WHERE protein_accession IS NOT NULL
            '''
        )
        return [x[0] for x in self.cursor.fetchall()]
    
    def get_organism_data_for_all_pubchem_uniprot_ids(self) -> Dict[str, List[str]]:
        pubchem_uniprot_ids = self.get_uniprot_ids_with_data() # get all unique uniprot IDs in the pubchem database
        uniprot_query = UniProtQuery() # open a connection to UniProt database
        out = uniprot_query.get_all_data_for_uniprot_ids(tuple(pubchem_uniprot_ids)) # get organism data for each \
            # uniprot ID
        uniprot_query.connection.close() # close the connection to the UniProt database
        return out
        
        
        

if __name__ == '__main__':
    pubchem_query = PubChemQuery()
    # pubchem_query.get_bioassays_from_uniprot_id('Q92952')
    # pubchem_query.get_bioassays_from_uniprot_id('Q6W5P4')
    # print(pubchem_query.get_represented_uniprot_ids())
    # print(pubchem_query.get_uniprot_ids_with_data())
    # print(len(pubchem_query.get_bioassay_protein_accessions()))
    pubchem_query.get_organism_data_for_all_pubchem_uniprot_ids()