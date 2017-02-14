from setuptools import setup
from sys import version_info, exit


if not (version_info.major == 3 and version_info.minor > 5):
    exit('The rxncon framework requires Python version 3.5 or higher. Please upgrade.')


setup(
    name='rxncon',
    packages=['rxncon'],
    version='2.0b1',
    description='The reaction-contingency framework for cellular signalling processes.',
    author='The rxncon group @ Humboldt University Berlin',
    author_email='jesper.romers@hu-berlin.de',
    url='https://github.com/rxncon/rxncon',
    download_url='https://github.com/rxncon/rxncon/tarball/2.0b1',
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
