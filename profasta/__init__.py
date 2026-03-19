"""ProFASTA - A Python library for working with protein containing FASTA files."""

from . import io, parser
from .db import DatabaseEntry, ProteinDatabase
from .decoy import create_decoy_db

__author__ = "David M. Hollenstein"
__license__ = "MIT"
__version__ = "0.0.5"
__all__ = [
    "DatabaseEntry",
    "ProteinDatabase",
    "create_decoy_db",
    "parser",
    "io",
]
