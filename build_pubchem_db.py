import json
import duckdb
import pyarrow as pa
import copy


def read_file(json_path:str) -> dict:
    with open(json_path) as f:
        return json.load(f)

def format_file(file_json:dict) -> pa.lib.Table:
    # Extract from file data
    data_copy = copy.deepcopy(file_json['PC_AssaySubmit']['data'])
    
    ########## Format data_copy ##########
    pylist = []
    for sid_entry in data_copy:
        sid_results = sid_entry.pop('data') # List of dicts

        # Extract useful data for each TID
        for i, tid_data in enumerate(sid_results):
            tid_data.update({str(tid_data.pop('tid')): # use string TID as dict key
                             list(tid_data.pop('value').values())[0]}) # list of one element
            sid_results[i] = tid_data # over-write

        # Convert list of dictionaries to single dict
        sid_results = {k: v for d in sid_results for k, v in d.items()}

        sid_entry.update(sid_results)
        pylist.append(sid_entry)
        
    data_table = pa.Table.from_pylist(pylist) # convert to table
        
    ########## Join data_copy w/ results_table on TIDs ##########
    exclude_names = ['sid', 'version', 'outcome', 'rank'] # non-TID column names
    # List of TID column names in data_table
    old_names = [int(x) for x in data_table.column_names if x not in exclude_names]
    # List of names mapped to TIDs in results_table
    results_table = pa.Table.from_pylist(file_json['PC_AssaySubmit']['assay']['descr']['results']) # Used for DuckDB
    new_names = [x[1] for x in duckdb.sql(f'SELECT * FROM results_table WHERE tid IN {tuple(old_names)}').fetchall()]
    
    return data_table.rename_columns(exclude_names + new_names)


if __name__ == '__main__':
    formatted_table = format_file(read_file('1865.json'))
    print(duckdb.sql('SELECT * FROM formatted_table').df())
