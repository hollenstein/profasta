"""Module for managing protein entries derived from FASTA files.

This module provides a ProteinDatabase class for managing protein entries derived from
FASTA files. The database can be used to import FASTA files, access the protein entries
by their identifier, and write the protein entries back to a FASTA file.

Classes:
    AbstractDatabaseEntry (Protocol): Interface for representing a protein entry.
    DatabaseEntry: Representation of a protein entry.
    ProteinDatabase: A database for managing protein entries derived from FASTA files.
"""

from __future__ import annotations
from dataclasses import dataclass
import logging
from pathlib import Path
from typing import Any, Iterator, Optional, Protocol

from profasta.parser import get_parser, get_writer
import profasta.io

logger = logging.getLogger(__name__)


class AbstractDatabaseEntry(Protocol):
    """A protein entry derived from a protein record in a FASTA file."""

    identifier: str
    header: str
    sequence: str
    header_fields: dict[str, str]


@dataclass
class DatabaseEntry:
    """A protein entry derived from a protein record in a FASTA file."""

    identifier: str
    header: str
    sequence: str
    header_fields: dict[str, str]


class ProteinDatabase:
    """A database for storing protein entries derived from FASTA files.

    This class provides functionality for managing a collection of protein entries
    derived from FASTA files. It allows for importing protein entries from FASTA files,
    adding new entries, and exporting entries back to FASTA format.

    The Python indexing operator `[]` provides access to the protein entries by their
    identifier. Iterating over the database yields the identifiers of all protein
    entries.

    Attributes:
        db: Dictionary mapping protein identifiers to protein entries.
        added_fasta_files: List of FASTA files that have been imported into the
            database.
        skipped_fasta_entries: Dictionary mapping added FASTA names to lists of FASTA
            entry headers that could not be parsed by the header parser, and thus were
            not added to the database.
    """

    db: dict[str, AbstractDatabaseEntry]
    added_fasta_files: list[str]
    skipped_fasta_entries: dict[str, list]

    def __init__(self):
        self.db = {}
        self.added_fasta_files = []
        self.skipped_fasta_entries = {}

    def add_fasta(
        self,
        path: str,
        header_parser: str,
        fasta_name: Optional[str] = None,
        overwrite: bool = False,
        skip_invalid: bool = False,
    ):
        """Add protein entries from a FASTA file to the database.

        Args:
            path: The path to the FASTA file.
            header_parser: The name of the parser to use for parsing the FASTA headers.
                The parser must be registered in the global parser registry.
            fasta_name: Optional name for the FASTA file. If not provided, the filename
                will be used instead.
            overwrite: If True, overwrite an existing entry with the same identifier.
                If False and an entry with the same identifier already exists, a
                KeyError will be raised.
            skip_invalid: If True, entries with a non-parsable header are skipped. If
                False, a ValueError is raised when an entry is encountered which header
                could not be parsed by the header_parser. Headers of skipped entries are
                stored in the skipped_fasta_entries attribute.
        """
        fasta_name = fasta_name if fasta_name is not None else Path(path).name
        parser = get_parser(header_parser)
        parsed_protein_entries: list[DatabaseEntry] = []
        skipped_entry_headers: list[str] = []

        with open(path, "r") as file:
            for fasta_record in profasta.io.parse_fasta(file):
                try:
                    parsed_header = parser.parse(fasta_record.header)
                except ValueError as error:
                    if skip_invalid:
                        skipped_entry_headers.append(fasta_record.header)
                        continue
                    else:
                        raise ValueError(
                            f"FASTA header could not be parsed with the "
                            f"'{header_parser}' parser: '{fasta_record.header}'"
                        ) from error
                protein_entry = DatabaseEntry(
                    parsed_header.identifier,
                    parsed_header.header,
                    fasta_record.sequence,
                    parsed_header.header_fields,
                )
                parsed_protein_entries.append(protein_entry)

        self.added_fasta_files.append(fasta_name)
        self.skipped_fasta_entries[fasta_name] = skipped_entry_headers
        for protein_entry in parsed_protein_entries:
            self.add_entry(protein_entry, overwrite)

        if skipped_entry_headers:
            num_skipped = len(skipped_entry_headers)
            num_total = num_skipped + len(parsed_protein_entries)
            logger.warning(
                f"Skipped {num_skipped}/{num_total} entries while adding "
                f"'{fasta_name}' to a ProteinDatabase because their headers could not "
                f"be parsed:"
            )

    def add_entry(self, protein_entry: AbstractDatabaseEntry, overwrite: bool = False):
        """Add a protein entry to the database.

        Args:
            protein_entry: The protein entry to add.
            overwrite: If True, overwrite an existing entry with the same identifier.
                If False and an entry with the same identifier already exists, a
                KeyError will be raised.
        """
        if protein_entry.identifier in self.db and not overwrite:
            raise KeyError(
                f"Identifier {protein_entry.identifier} already in database."
            )

        self.db[protein_entry.identifier] = protein_entry

    def write_fasta(
        self,
        path,
        append: bool = False,
        header_writer: Optional[str] = None,
        line_width: int = 60,
    ):
        """Write all protein entries in the database to a FASTA file.

        Args:
            path: The path to write the FASTA file to.
            append: If False, the file is created or overwritten. If True, the entries
                are appended to an existing file. The default value is False.
            header_writer: The name of the writer to use for generating the FASTA header
                strings from the database entries. If None, the original header strings
                are written to the FASTA file.
            line_width: The number of sequence characters per line, the default value is
                60. If -1, the sequence is not split into multiple lines.
        """
        if header_writer is None:
            fasta_records = list(self.db.values())
        else:
            fasta_records = []
            writer = get_writer(header_writer)
            for protein_entry in self.db.values():
                header = writer.write(protein_entry)

                fasta_records.append(
                    DatabaseEntry("", header, protein_entry.sequence, {})
                )
        file_open_mode = "a" if append else "w"
        with open(path, file_open_mode) as file:
            profasta.io.write_fasta(file, fasta_records, line_width)

    def get(self, identifier: str, default: Any = None) -> DatabaseEntry | Any:
        """Get a protein entry by its identifier or return a default value."""
        return self.db.get(identifier, default)

    def keys(self):
        return self.db.keys()

    def values(self):
        return self.db.values()

    def items(self):
        return self.db.items()

    def __getitem__(self, identifier) -> AbstractDatabaseEntry:
        return self.db[identifier]

    def __contains__(self, identifier) -> bool:
        return identifier in self.db

    def __iter__(self) -> Iterator[str]:
        return iter(self.db)
