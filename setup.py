from setuptools import setup, find_packages
from mag_annotator import __version__ as version

__author__ = 'shafferm'
__version__ = version

setup(
      name="MAGotator",
      version=__version__,
      scripts=['scripts/MAGotator.py', 'scripts/annotate_genes.py', 'scripts/prepare_databases.py',
               'scripts/make_genome_summary.py', 'scripts/set_database_locations.py'],
      packages=find_packages(),
      package_data={'mag_annotator': ['DATABASE_LOCATIONS']},
      install_requires=['scikit-bio', 'pandas', 'networkx'],
      description="Annotate contigs/bins from metagenomic assemblies and create predicted metabolisms",
      author="Michael Shaffer",
      author_email='michael.t.shaffer@colostate.edu',
      url="https://github.com/shafferm/MAGotator/",
      download_url="https://github.com/shafferm/MAGotator/tarball/%s" % __version__
)
