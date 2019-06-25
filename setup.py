from setuptools import setup, find_packages
from checkMetab import __version__ as version

__author__ = 'shafferm'
__version__ = version

setup(
      name="checkMetab",
      version=__version__,
      scripts=['scripts/annotate_genes.py', 'scripts/prepare_databases.py', 'scripts/make_genome_summary.py',
               'scripts/set_database_locations.py'],
      packages=find_packages(),
      package_data={'checkMetab': ['DATABASE_LOCATIONS']},
      install_requires=['scikit-bio', 'pandas', 'networkx'],
      description="Annotate contigs/bins from metagenomic assemblies and create predicted metabolisms",
      author="Michael Shaffer",
      author_email='michael.t.shaffer@colostate.edu',
      url="https://github.com/shafferm/checkMetab/",
      download_url="https://github.com/shafferm/checkMetab/tarball/%s" % __version__
)
