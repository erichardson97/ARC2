import requests
import subprocess

def download_igblast():
    out = requests.get('https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.22.0-x64-linux.tar.gz')
    with open('igblast.linux.tar.gz', 'wb') as k:
        k.write(out.content)
    subprocess.call(f'tar -xzvf igblast.linux.tar.gz', shell=True)