import os
import shutil
from ftplib import FTP
from functools import wraps
from tqdm import tqdm


class PubChemFTP():
    def __init__(self, out_dir:str, overwrite:bool=False):
        self.__out_dir = out_dir
        self.__overwrite = overwrite
        
        self._make_dir()
        
        # FTP connection details
        ftp_host = 'ftp.ncbi.nlm.nih.gov'
        self.__ftp_user = 'anonymous'
        self.__ftp_password = ''
        
        # Connect to the FTP server
        self.__ftp = FTP(ftp_host)
        self.__ftp.login(self.__ftp_user, self.__ftp_password)
    
    def download_all(self, verbose:bool=True):
        # Download protein target data 'protein2xrefs.gz'
        if verbose:
            print('Downloading protein target data...')
        try:
            self.download_protein_target_data(verbose=verbose)
        except Exception as e:
            print(f'Error downloading protein target data: {e}')
        
        # Download substance SDFs
        if verbose:
            print('Downloading substance SDFs...')
        try:
            self.download_substance_sdfs(verbose=verbose)
        except Exception as e:
            print(f'Error downloading substance SDFs: {e}')
        
        # Download bioassay JSONs
        if verbose:
            print('Downloading bioassay JSONs...')
        try:
            self.download_bioassay_jsons(verbose=verbose)
        except Exception as e:
            print(f'Error downloading bioassay JSONs: {e}')
            
    def download_protein_target_data(self, verbose:bool=True) -> None:
        filename = 'protein2xrefs.gz'

        # Change to the target directory
        ftp_directory = 'pubchem/Target/'
        self.__ftp.cwd('/')
        self.__ftp.cwd(ftp_directory)
        
        protein_target_out_dir = os.path.join(os.getcwd(), self.__out_dir, ftp_directory)
        os.makedirs(protein_target_out_dir)
        file_path = os.path.join(protein_target_out_dir, filename)
        
        with open(file_path, 'wb') as file:
            self.__ftp.retrbinary('RETR ' + filename, file.write)
        if verbose:
            print(f'Downloaded: {filename}')

    def download_substance_sdfs(self, verbose:bool=True) -> None:
        # FTP directory and file details
        ftp_directory = 'pubchem/Substance/CURRENT-Full/SDF/'
        self.__ftp.cwd('/')
        self.__ftp.cwd(ftp_directory)
        
        substance_sdf_out_dir = os.path.join(os.getcwd(), self.__out_dir, ftp_directory)
        os.makedirs(substance_sdf_out_dir)

        # Retrieve a list of all file names in the directory
        filenames = [file_name for file_name in self.__ftp.nlst() if not file_name.endswith('.sdf.gz.md5')] # exclude \
            # MD5 checksum files, but retain the README files

        # Download each file
        for filename in (tqdm(filenames) if verbose else filenames):
            try:
                file_path = os.path.join(substance_sdf_out_dir, filename)
                with open(file_path, 'wb') as file:
                    self.__ftp.retrbinary('RETR ' + filename, file.write)
                if verbose:
                    print(f"Downloaded: {filename}")
            except Exception as e:
                print(f'Error downloading {filename}: {e}')
        
    def download_bioassay_jsons(self, verbose:bool=True) -> None:
        # Change to the target directory
        ftp_directory = 'pubchem/Bioassay/JSON/'
        self.__ftp.cwd('/')
        self.__ftp.cwd(ftp_directory)
        
        bioassay_json_out_dir = os.path.join(os.getcwd(), self.__out_dir, ftp_directory)
        os.makedirs(bioassay_json_out_dir)

        # Retrieve a list of all file names in the directory
        filenames = self.__ftp.nlst()

        # Download each file
        for filename in (tqdm(filenames) if verbose else filenames):
            try:
                file_path = os.path.join(bioassay_json_out_dir, filename)
                with open(file_path, 'wb') as file:
                    self.__ftp.retrbinary('RETR ' + filename, file.write)
                if verbose:
                    print(f'Downloaded: {filename}')
            except Exception as e:
                print(f'Error downloading {filename}: {e}')
                
    def _make_dir(self):
        """
        Conditionally create a directory for the output files.

        :raises ValueError: Overwrite is False but directory is not empty
        """
        if not os.path.exists(self.__out_dir): # create directory if it doesn't exist
            os.makedirs(self.__out_dir)
        elif len(os.listdir(self.__out_dir)) > 0 and self.__overwrite == False: # raise error if directory exists and \
            # overwrite is False
            raise ValueError(f'Directory {self.__out_dir} is not empty. Set overwrite=True to overwrite files. Note \
                that this will delete all files in {self.__out_dir}.')
        else: # delete directory and recreate if overwrite is True
            shutil.rmtree(self.__out_dir)
            os.makedirs(self.__out_dir)
        

if __name__ == '__main__':
    pc_ftp = PubChemFTP('pubchem_ftp_data', overwrite=True)
    pc_ftp.download_all(verbose=True)
