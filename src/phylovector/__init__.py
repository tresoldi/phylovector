"""
phylovector __init__.py
"""

# Version of the phylovector package
__version__ = "0.1"
__author__ = "Tiago Tresoldi"
__email__ = "tiago.tresoldi@lingfil.uu.se"

# Import from local modules
from .phylovector import tree2vector, vector2tree
from .newick import _parse, sorted_newick

# Build the namespace
__all__ = [
    "tree2vector",
    "vector2tree",
    "_parse", # TODO: needed in the common space?
    "sorted_newick",
]
