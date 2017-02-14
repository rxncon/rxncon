from setuptools import setup, find_packages
from sys import version_info, exit
from platform import platform


def install_pyeda_from_alternate_source():
    reply = input('This installer will install a binary version of the PyEDA package\nsince you run on Windows.'
                  ' Please confirm by typing \'yes\':\n')

    if reply != 'yes':
        exit(1)

    # find out 32-bit or 64-bit
    # find out python 3.5 or python 3.6
    # construct pyeda.xxxxx.whl filename
    # download file (curl, wget, anything)
    # put into temporary folder
    # shell out and do 'pip install pyeda-xxx.whl" ? dunno which rights required pip install --user??

    pass


if not (version_info.major == 3 and version_info.minor >= 5):
    exit('The rxncon framework requires Python version 3.5 or higher. Please upgrade.')

windows = platform().split('-')[0] == 'Windows'

if windows:
    install_pyeda_from_alternate_source()

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
    version='2.0b4',
    description='The reaction-contingency framework for cellular signalling processes.',
    author='The rxncon group @ Humboldt University Berlin',
    author_email='jesper.romers@hu-berlin.de',
    url='https://github.com/rxncon/rxncon',
    download_url='https://github.com/rxncon/rxncon/tarball/2.0b4',
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
