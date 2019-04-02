from setuptools import setup, find_packages
from checkMetab import __version__ as version

__author__ = 'shafferm'
__version__ = version

setup(
      name="checkMetab",
      version=__version__,
      scripts=['scripts/annotate_genes.py', 'scripts/generate_database_form.py', 'scripts/check_metab.py'],
      packages=find_packages(),
      install_requires=['scikit-bio', 'pandas', 'networkx'],
      description="Annotate contigs/bins from metagenomic assemblies and create predicted metabolisms",
      author="Michael Shaffer",
      author_email='michael.t.shaffer@colostate.edu',
      url="https://github.com/shafferm/checkMetab/",
      download_url="https://github.com/shafferm/checkMetab/tarball/%s" % __version__
)
