from autochem import PubChem
from ChemiCure.chemicure import ChemiCure
import pandas as pd

target_id_1 = 'P50129' # chembl_id = CHEMBL2490
target_id_2 = 'A0AVT1'
pubchem = PubChem([target_id_1, target_id_2])
# pubchem = PubChem([target_id_2])

# df = pubchem.df

# ChemiCure()
# df.to_csv('test_chembl.csv', index=False)
# print(df.to_sql())