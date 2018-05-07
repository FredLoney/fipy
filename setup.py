import os
import re
from setuptools import (setup, find_packages)


class InstallError(Exception):
    """reactome fipy installation error."""
    pass


def version(package):
    """
    :return: the package version as listed in the package `__init.py__`
        `__version__` variable.
    """
    # The version string parser.
    REGEXP = re.compile("""
       __version__   # The version variable
       \s*=\s*       # Assignment
       ['\"]         # Leading quote
       (.+)          # The version string capture group
       ['\"]         # Trailing quote
    """, re.VERBOSE)

    with open(os.path.join(package, '__init__.py')) as f:
       match = REGEXP.search(f.read())
       if not match:
           raise InstallError("The reactome fipy __version__ variable"
                              " was not found")
       return match.group(1)


def requires():
    with open('requirements.txt') as f:
        return f.read().splitlines()


def readme():
    with open("README.rst") as f:
        return f.read()


setup(
    name = 'reactome-fipy',
    version = version('reactome/fipy'),
    author = 'Oregon Health & Science University',
    author_email = 'loneyf@ohsu.edu',
    platforms = 'Any',
    license = 'MIT',
    keywords = 'Reactome Cytoscape pathway enrichment',
    packages = find_packages(exclude=['test**']),
    url = 'http://reactome-fipy.readthedocs.org/en/latest/',
    description = 'Rectome FI CyREST Facade',
    long_description = readme(),
    classifiers = [
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
    ],
    install_requires = requires()
)
