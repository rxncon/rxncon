from setuptools import setup, find_packages
from sys import version_info, exit


if not (version_info.major == 3 and version_info.minor > 5):
    exit('The rxncon framework requires Python version 3.5 or higher. Please upgrade.')


setup(
    name='rxncon',
    packages=find_packages(),
    package_data={
        '': ['*.xls', '*.tsv', '*.xgmml']
    },
    scripts=[
        'rxncon2bngl.py',
        'rxncon2boolnet.py',
        'rxncon2reactiongraph.py',
        'rxncon2regulatorygraph.py'
    ],
    version='2.0b3',
    description='The reaction-contingency framework for cellular signalling processes.',
    author='The rxncon group @ Humboldt University Berlin',
    author_email='jesper.romers@hu-berlin.de',
    url='https://github.com/rxncon/rxncon',
    download_url='https://github.com/rxncon/rxncon/tarball/2.0b3',
    license='LGPL',
    keywords=['sysbio', 'signalling', 'systems biology'],
    classifiers=[],
    install_requires=[
        'pytest',
        'numpy',
        'scipy',
        'pyeda',
        'click',
        'click_log',
        'colorama',
        'xlrd',
        'networkx',

    ]
)
