import setuptools
#from ARC.download_data import *
from setuptools.command.install import install
import shutil
import os

class PostInstallCommand(install):
    def run(self):
        install.run(self)
      #  print(f'Downloading and processing relevant IMGT and MRO data.')
      #  data = DataDownloader()
      #  data.download_MRO_tsv()
      #  data.download_IG_TR_databases()
      #  for file in os.listdir(os.path.join(data.package_directory, 'data')):
      #      if os.path.exists(os.path.join(self.install_lib, 'ARC', 'data', file)):
      #          if os.path.isdir(os.path.join(self.install_lib, 'ARC', 'data', file)):
      #              shutil.rmtree(os.path.join(self.install_lib, 'ARC', 'data', file))
      #              shutil.copytree(os.path.join(data.package_directory, 'data', file), os.path.join(self.install_lib, 'ARC', 'data', file))
      #          else:
      #              os.remove(os.path.join(self.install_lib, 'ARC', 'data', file))
      #              shutil.copyfile(os.path.join(data.package_directory, 'data', file), os.path.join(self.install_lib, 'ARC', 'data', file))
      #      
      #      elif os.path.isdir(os.path.join(data.package_directory, 'data', file)):
      #          shutil.copytree(os.path.join(data.package_directory, 'data', file), os.path.join(self.install_lib, 'ARC', 'data', file))
      #      else:
      #          shutil.copyfile(os.path.join(data.package_directory, 'data', file), os.path.join(self.install_lib, 'ARC', 'data', file))
        

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bio-arc",
    version="0.2",
    author="Austin Crinklaw",
    author_email="erichardson@lji.org",
    description="Antigen Receptor Classifier (II)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/erichardson97/ARC2",
    packages=setuptools.find_packages(),
    package_data={
        'ARC': [
            'data/*', 'data/blastdb/*', 'data/HMMs/*', 'data/IgNAR/**',
            'data/MHC_HMMs/*', 'tests/', 'muscle', 'data/imgt/**', 'data/imgt/blast_fasta/*'
        ]
    },
    include_package_data=True,
    install_requires=[
        'pandas',
        'biopython',
        'numpy',
        'tqdm',
        'requests',
        'beautifulsoup4'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
    ])
