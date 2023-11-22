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
            SELECT assay_data, smiles, bioassay_id
            FROM (
                SELECT bioassay_id, json_array_elements(assay_data) AS assay_data, json_array_elements(assay_data)->>'sid' AS sid
                FROM bioassay
                JOIN temp_assay_ids ON bioassay_id = id
            ) AS extracted_sids
            JOIN substance ON sid::integer = substance_id
        ''')
        data = cursor.fetchall()
        cursor.execute('DROP TABLE temp_assay_ids')
        
        df = pd.DataFrame([x[0] for x in data])  # assay data
        df.insert(0, 'smiles', [x[1] for x in data])
        df.insert(0, 'bioassay_id', [x[2] for x in data])
        df.insert(0, 'uniprot_id', uniprot_id)
        dfs.append(df)
        
    cursor.close()
    connection.close()

    dfs[0].to_csv('ache_assay_data.csv', index=False)
    dfs[1].to_csv('bche_assay_data.csv', index=False)

