from chemcurate import Chembl
from ChemiCure.chemicure import ChemiCure
import pandas as pd

target_id_1 = 'P50129' # chembl_id = CHEMBL2490
target_id_2 = 'Q13936'
chembl = Chembl([target_id_1, target_id_2])

df = chembl.df

# ChemiCure()
# df.to_csv('test_chembl.csv', index=False)
# print(df.to_sql())