"""
setup.py for the `phylovector` library

Requirements are listed in `requirements.txt`.
"""


# Import Python standard libraries
from setuptools import setup, find_packages
import pathlib

# The directory containing this file
LOCAL_PATH = pathlib.Path(__file__).parent

# The text of the README file
README_FILE = (LOCAL_PATH / "README.md").read_text(encoding="utf-8")

# Load requirements, so they are listed in a single place
with open("requirements.txt", encoding="utf-8") as fp:
    install_requires = [dep.strip() for dep in fp.readlines()]

# This call to setup() does all the work
setup(
    author="Tiago Tresoldi",
    author_email="tiago.tresoldi@lingfil.uu.se",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Software Development :: Libraries",
    ],
    description="Library for the vector representation of phylogenetic trees",
    entry_points={"console_scripts": ["phylovector=phylovector.__main__:main"]},
    extras_require={
        "paper" : ["ngesh"],
        "dev": ["black", "flake8", "twine", "wheel"],
        "test": ["pytest"],
    },
    include_package_data=True,
    install_requires=install_requires,
    keywords=["phylogenetic tree", "phylogenetics", "optimization"],
    license="MIT",
    long_description=README_FILE,
    long_description_content_type="text/markdown",
    name="phylovector",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.6",
    test_suite="tests",
    tests_require=[],
    url="https://github.com/tresoldi/phylovector",
    version="0.1",  # remember to sync with __init__.py
    zip_safe=False,
)
