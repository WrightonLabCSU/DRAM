from setuptools import setup, find_packages
from mag_annotator import __version__ as version

__author__ = 'shafferm'
__version__ = version

setup(
    name="DRAM",
    version=__version__,
    scripts=['scripts/DRAM.py'],
    packages=find_packages(),
    package_data={'mag_annotator': ['CONFIG']},
    install_requires=['scikit-bio', 'pandas', 'altair', 'sqlalchemy', 'networkx', 'openpyxl', 'numpy'],
    description="Annotate contigs/bins from metagenomic assemblies and create predicted metabolisms",
    author="Michael Shaffer",
    author_email='michael.t.shaffer@colostate.edu',
    url="https://github.com/shafferm/DRAM/",
    download_url="https://github.com/shafferm/DRAM/tarball/%s" % __version__
)
