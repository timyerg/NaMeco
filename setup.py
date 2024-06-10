from setuptools import setup

setup(
    name='nameco',
    version='0.1.1',    
    description='Pipeline for the Nanopore 16S long read clustering and taxonomy classification',
    url='https://github.com/timyerg/NaMeco',
    author='Timur Yergaliyev',
    author_email='timyerg@gmail.com',
    license='Apache-2.0',
    packages=['nameco'],
    install_requires=['pandas>=2.2',
		      'python>=3.10.8',
		      'chopper>=0.7',
		      'racon>=1.5',
		      'medaka>=1.11',
		      'minimap2>=2.28',
		      'biopython>=1.83',
		      'matplotlib>=3.8',
		      'blast>=2.15',
		      'spoa>=4.1',
		      'scikit-learn>=1.5',
		      'umap-learn>=0.5'],)
