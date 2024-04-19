"""ProFASTA - A Python library for working with protein containing FASTA files."""

from profasta.db import DatabaseEntry, ProteinDatabase
from profasta.decoy import create_decoy_db
import profasta.parser
import profasta.io


__author__ = "David M. Hollenstein"
__license__ = "MIT"
__version__ = "0.0.5"
