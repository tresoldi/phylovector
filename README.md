# phylovector

Library for the vector representation of phylogenetic trees

`philovector` implements the method for representing phylogenetic trees as vectors (i.e., one-dimensional arrays) as described in Tresoldi (forth.). The library was developed as a high-level Python script, with code as close as possible to pseudo-code, in order to facilitate the conversion to other programming languages. Note that, while offered as a normal Python package, the code is designed to be easily reused in other projects without installing 3rd-party libraries; as such, it is entirely contained in a single script that can be copied to other projects and has no external dependecies with the exception of the `ete3` and `networkx` libraries (and its dependencies).

## Installation

To execute the demonstrations and experiments of the paper, the `extra_requirements` listed in `setup.py` must be installed as well. This is usually done with `pip install -e .[paper]` (note that you might need to escape the brackets depending on you shell configuration).