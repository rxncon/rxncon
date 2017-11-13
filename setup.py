
from setuptools import setup, find_packages
from sys import version_info, exit
from platform import platform, architecture
from pip import main as pip_main
from urllib.request import urlopen
from urllib.error import URLError
from ssl import CERT_NONE, create_default_context


def install_precompiled_windows_pyeda():
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
        system_version = architecture()[0]
        major, minor, micro, releaselevel, serial = version_info

        pyeda_str = 'pyeda-0.28.0-cp{0}{1}-cp{0}{1}m-win'.format(major, minor)

        if system_version == "32bit":
            pyeda_str += '32.whl'
        elif system_version == '64bit':
            pyeda_str += '_amd64.whl'
        else:
            raise OSError('System version {} not known.'.format(system_version))

        source_path = 'http://rumo.biologie.hu-berlin.de/rxncon_downloads/'
        try:
            urlopen(source_path + pyeda_str, context=ssl_context())
            pip_main(['install', source_path + pyeda_str])
        except URLError:
            exit(1)


if not (version_info.major == 3 and version_info.minor >= 5):
    exit('The rxncon framework requires Python version 3.5 or higher. Please upgrade.')


windows = platform().split('-')[0] == 'Windows'


if windows:
    install_precompiled_windows_pyeda()


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
        'rxncon2regulatorygraph.py',
        'rxncon2srgraph.py'
    ],
    version='2.0b17',
    description='The reaction-contingency framework for cellular signalling processes.',
    author='The rxncon group @ Humboldt University Berlin',
    author_email='jesper.romers@hu-berlin.de',
    url='https://github.com/rxncon/rxncon',
    license='LGPL',
    keywords=['sysbio', 'signalling', 'systems biology'],
    classifiers=[],
    install_requires=[
        'pytest',
        'numpy',
        'scipy',
        'pyeda',
        'click',
        'click_log==0.1.8',
        'colorama',
        'xlrd',
        'networkx',
    ]
)
