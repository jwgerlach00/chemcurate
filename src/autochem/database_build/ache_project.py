import requests
import psycopg2
import pandas as pd


def get_ache_bche_assay_ids() -> dict:
    assay_ids = {}
    for uniprot_id in ['P22303', 'P06276']:
        _url_stem = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
        url = f'{_url_stem}/bioassay/target/ProteinName/{uniprot_id}/aids/JSON'
        
        json_response = requests.get(url).json()
        if 'IdentifierList' not in json_response:
            continue
        
        assay_ids[uniprot_id] = [(x,) for x in requests.get(url).json()['IdentifierList']['AID']]
    return assay_ids


if __name__ == '__main__':
    assay_ids = get_ache_bche_assay_ids()
    
    ########## Connect to DB ##########
    host = 'localhost'
    database = 'pubchem'
    user = 'postgres'
    password = 'Coll@bor@tions2020'
    # port = '5432'
    connection = psycopg2.connect(host=host, database=database, user=user, password=password)
    cursor = connection.cursor()
    
    dfs = []
    for uniprot_id in assay_ids:
        cursor.execute('CREATE TEMPORARY TABLE temp_assay_ids (id INT)')

        # Insert data into the temporary table
        cursor.executemany('INSERT INTO temp_assay_ids VALUES (%s)', assay_ids[uniprot_id])  # only look at ache for now
        connection.commit()
        
        cursor.execute('SELECT * FROM temp_assay_ids')
        cursor.execute('''
            SELECT assay_data, smiles
            FROM (
                SELECT *
                FROM (
                    SELECT json_array_elements(assay_data) AS assay_data, json_array_elements(assay_data)->>'sid' AS sid
                    FROM bioassay
                    JOIN temp_assay_ids ON bioassay_id = id
                ) AS extracted_sids
                JOIN substance ON sid::integer = substance_id
            ) AS test;
        ''')
        data = cursor.fetchall()
        cursor.execute('DROP TABLE temp_assay_ids')
        
        assay_data = [x[0] for x in data]
        smiles = [x[1] for x in data]
        df = pd.DataFrame(assay_data)
        df.insert(0, 'smiles', smiles)
        dfs.append(df)
        
    cursor.close()
    connection.close()

    dfs[0].to_csv('ache_assay_data.csv', index=False)
    dfs[1].to_csv('bche_assay_data.csv', index=False)

