"""
phylovector __init__.py
"""

# Version of the phylovector package
__version__ = "0.1"
__author__ = "Tiago Tresoldi"
__email__ = "tiago.tresoldi@lingfil.uu.se"

# Import from local modules
from .phylovector import tree2vector, vector2tree
from .common import sorted_newick, start_vector
from .likelihood import method1, primate

def vector_length(num_taxa):
    return int((pow(num_taxa, 2) + num_taxa) / 2.0)

# Build the namespace
__all__ = [
    "sorted_newick",
    "tree2vector",
    "vector2tree",
    "method1",
    "start_vector",
    "primate",
    "vector_length",
]
