
from setuptools import setup, find_packages
from sys import version_info, exit
from platform import platform, architecture
from pip import main as pip_main
from urllib.request import urlopen
from urllib.error import URLError
from ssl import CERT_NONE, create_default_context


def install_pyeda_from_alternate_source():
    def ssl_context():
        """
        create an ssl context, where the certificate is not verified

        Notes:
            This is necessary because http://rumo.biologie.hu-berlin.de has no valid certificate.
        Returns:
            SSLContext

        """
        context = create_default_context()
        context.check_hostname = False
        context.verify_mode = CERT_NONE
        return context

    try:
        import pyeda
    except ImportError:
        reply = input('This installer will install a binary version of the PyEDA package\nsince you run on Windows.'
                      ' Please confirm by typing \'yes\':\n')

        if reply != 'yes':
            exit(1)

        # find out 32-bit or 64-bit
        system_version = architecture()[0]
        print('System Version: {}'.format(system_version))
        # find out python 3.5 or python 3.6
        major, minor, micro, releaselevel, serial = version_info
        print('Python Version: {0}.{1}.{2}'.format(major, minor, micro))
        if major < 3 or minor < 5:
            raise "Need python 3.5 or greater"

        # construct pyeda.xxxxx.whl filename
        pyeda_str = 'pyeda-0.28.0-cp{0}{1}-cp{0}{1}m-win'.format(major, minor)
        if system_version == "32bit":
            pyeda_str += '32.whl'
        elif system_version == '64bit':
            pyeda_str += '_amd64.whl'
        print('PyEDA: {}'.format(pyeda_str))
        # install pyeda-xxx.whl with pip
        SOURCE_PATH = 'http://rumo.biologie.hu-berlin.de/rxncon_downloads/'
        try:
            urlopen(SOURCE_PATH + pyeda_str, context=ssl_context())
            pip_main(['install', SOURCE_PATH + pyeda_str])
        except URLError:
            print('No pre-compiled PyEDA version for your system configuration Windows architecture: {0} Python version: {1} available.'.format(system_version, ".".join([str(major), str(minor), str(micro)])))


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

