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
    entry_points={'console_scripts': ['nameco=run_nameco:run_pipeline',],},
     )
