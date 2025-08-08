# ARC2
ARC with added species assignment either by  HMMER (significantly less accurate; default) or BLAST (slow). Add the -s flag for BLAST.
[July 2025] - Added a separate more efficient classifier file operating on the FASTA file - use -m fasta.
            - Added mmseqs functionality - use --speedy

## Requirements

Linux OS
HMMER3
Python 3+
Python packages: Pandas, BioPython

```shell
conda install -c conda-forge biopython -y
conda install -c bioconda hmmer=3.3.2 -y
```
