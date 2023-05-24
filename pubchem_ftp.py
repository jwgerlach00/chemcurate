import os
import shutil
from ftplib import FTP


def ftp_download(out_dir:str, overwrite:bool=False) -> None:
     
    
    if not os.path.exists(out_dir): # create directory if it doesn't exist
        os.makedirs(out_dir)
    elif len(os.listdir(out_dir)) > 0 and overwrite == False: # raise error if directory exists and overwrite is False
        raise ValueError(f'Directory {out_dir} is not empty. Set overwrite=True to overwrite files. \
            Note that this will delete all files in {out_dir}.')
    else: # delete directory and recreate if overwrite is True
        shutil.rmtree(out_dir)
        os.makedirs(out_dir)
        
    # FTP connection details
    ftp_host = 'ftp.ncbi.nlm.nih.gov'
    ftp_user = 'anonymous'
    ftp_password = ''

    # FTP directory and file details
    ftp_directory = '/pubchem/Substance/CURRENT-Full/SDF/'

    # Connect to the FTP server
    ftp = FTP(ftp_host)
    ftp.login(ftp_user, ftp_password)

    # Change to the target directory
    ftp.cwd(ftp_directory)

    # Retrieve a list of all file names in the directory
    file_names = [file_name for file_name in ftp.nlst() if not file_name.endswith('.sdf.gz.md5')] # exclude MD5 checksum
        # files, but retain the README files

    # Download each file
    for file_name in file_names:
        file_path = os.path.join(os.getcwd(), out_dir, file_name)
        with open(file_path, 'wb') as file:
            ftp.retrbinary('RETR ' + file_name, file.write)
        print(f"Downloaded: {file_name}")

    # Close the FTP connection
    ftp.quit()

if __name__ == '__main__':
    ftp_download('test_pubchem_ftp_download', overwrite=True)
