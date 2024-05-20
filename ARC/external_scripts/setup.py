import requests
import subprocess

def download_igblast():
    requests.get('https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.22.0-x64-linux.tar.gz')
    subprocess.call(f'tar -xzvf ncbi-igblast-1.22.0-x64-linux.tar.gz')