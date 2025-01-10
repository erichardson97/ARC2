import setuptools
from download_data import *

class PostInstallCommand(setuptools.command.install):
    def run(self):
        install.run(self)
        print(f'Downloading and processing relevant IMGT and MRO data.')
        data = DataDownloader()
        data.download_MRO_tsv()
        data.download_IG_TR_databases()


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
            'data/*', 'data/blastdb', 'data/HMMs', 'data/IgNAR',
            'data/MHC_HMMs', 'tests/*'
        ]
    },
    include_package_data=True,
    install_requires=[
        'pandas',
        'biopython',
        'numpy',
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
    ],
    cmdclass={"install": PostInstallCommand, },
    )
