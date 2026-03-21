"""ProFASTA - A Python library for working with protein containing FASTA files."""

from . import db, decoy, io, parser, validation
from .db import DatabaseEntry, ProteinDatabase
from .decoy import create_decoy_db, write_decoy_fasta

__author__ = "David M. Hollenstein"
__license__ = "MIT"
__version__ = "0.1.0"
__all__ = [
    "DatabaseEntry",
    "ProteinDatabase",
    "create_decoy_db",
    "write_decoy_fasta",
    # submodules
    "db",
    "decoy",
    "parser",
    "io",
    "validation",
]
