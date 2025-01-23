# ARC2
ARC with added species assignment either by BLAST (default; slow) or HMMER (significantly less accurate; fast).

## Requirements

Linux OS
HMMER3
NCBI IgBLAST
NCBI Blast+
Python 3+
Python packages: Pandas, BioPython, tqdm, bs4

```shell
conda install -c conda-forge biopython -y
conda install -c bioconda hmmer=3.3.2 -y
wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.22.0-x64-linux.tar.gz
tar -xvzf ncbi-igblast-1.22.0-x64-linux.tar.gz
## add to your .bashrc, move to existing bin path, or temporarily:
export PATH=$PATH:${PWD}/ncbi-igblast-1.22.0/bin
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
tar -xzvf ncbi-blast-2.16.0+-x64-linux.tar.gz
## add to your .bashrc, move to existing bin path, or temporarily:
export PATH=$PATH:${PWD}/ncbi-blast-2.16.0+/bin
```
