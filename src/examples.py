import os
from autochem.database_build import PubChemFTP, PubChemDB


FTP_ABSOLUTE_PATH = '/pubchem_ftp_data'


if __name__ == '__main__':
    ########## Download data from PubChem using the FTP server ##########
    pubchem_ftp = PubChemFTP(os.path.join(os.getcwd(), 'pubchem_ftp_data'), overwrite=True)
    
    # Download protein target xref, substance sdf, and bioassay json files
    pubchem_ftp.download_all(verbose=True)
    
    ########## Define paths to data downloaded from PubChem FTP site using autochem.database_build.PubChemFTP ##########
    pc_db = PubChemDB(
        f'{FTP_ABSOLUTE_PATH}/pubchem/Bioassay/JSON',
        f'{FTP_ABSOLUTE_PATH}/pubchem/Substance/CURRENT-Full/SDF',
        f'{FTP_ABSOLUTE_PATH}/pubchem/Target/protein2xrefs'
    )
    
    ########## Build the database using PostgreSQL ##########
    pass
    
    ########## Populate entries into the database ##########
    pc_db.repopulate_protein_target_table()
    pc_db.repopulate_substance_table()
    pc_db.repopulate_bioassay_table()
    # Container method for all of the above
    pc_db.build()
