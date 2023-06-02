from autochem.database_build import PubChemDB


if __name__ == '__main__':
    pc_db = PubChemDB(
        '/Users/collabpharma/Desktop/JSON',
        '/Users/collabpharma/Desktop/SDF',
        '/Users/collabpharma/Desktop/protein2xrefs'
    )
    pc_db.repopulate_bioassay_table()
