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

import logging
import pathlib
from dataclasses import dataclass
from os import PathLike
from typing import (
    Any,
    Callable,
    ItemsView,
    Iterator,
    KeysView,
    Optional,
    Protocol,
    ValuesView,
)

import profasta.io
from profasta.parser import get_parser, get_writer

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

    The class implements a dict-like interface:
    - Access entries by identifier: `db[identifier]`
    - Safe access with default: `db.get(identifier, default)`
    - Check membership: `identifier in db`
    - Iterate over identifiers: `for identifier in db`
    - Iterate over entries: `db.keys()`, `db.values()`, `db.items()`
    - Get database size: `len(db)`
    - Check if non-empty: `bool(db)` or `if db:`

    Attributes:
        added_fasta_files: List of FASTA files that have been imported into the
            database.
        skipped_fasta_entries: Dictionary mapping added FASTA names to lists of FASTA
            entry headers that could not be parsed by the header parser, and thus were
            not added to the database.
    """

    added_fasta_files: list[str]
    skipped_fasta_entries: dict[str, list[str]]

    def __init__(self):
        self._db = {}
        self.added_fasta_files = []
        self.skipped_fasta_entries = {}

    @classmethod
    def from_fasta(
        cls,
        *paths: PathLike,
        header_parser: str,
        overwrite: bool = False,
        skip_invalid: bool = False,
    ) -> "ProteinDatabase":
        """Create a ProteinDatabase from one or more FASTA files.

        A convenience constructor that creates a new database and immediately populates
        it by calling `.add_fasta` for each supplied path. All paths are parsed with the
        same `header_parser`.

        Args:
            *paths: One or more paths to FASTA files to import.
            header_parser: The name of the parser to use for parsing the FASTA headers.
                The parser must be registered in the global parser registry.
            overwrite: If True, overwrite existing entries with duplicate identifiers.
                If False and any entry from the file has an identifier already present
                in the database, a KeyError is raised before any entries are added.
            skip_invalid: If True, entries with a non-parsable header are skipped. If
                False, a ValueError is raised when an entry is encountered which header
                could not be parsed by the header_parser. Headers of skipped entries are
                stored in the skipped_fasta_entries attribute.

        Returns:
            A `ProteinDatabase` populated with all entries from the supplied FASTA files.

        Raises:
            KeyError: If `overwrite` is `False` and a duplicate identifier is found
                across the imported files.
            ValueError: If `skip_invalid` is `False` and a FASTA header cannot be parsed
                by `header_parser`.
        """
        db = cls()
        for path in paths:
            db.add_fasta(
                path,
                header_parser,
                overwrite=overwrite,
                skip_invalid=skip_invalid,
            )
        return db

    def add_fasta(
        self,
        path: PathLike,
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
            overwrite: If True, overwrite existing entries with duplicate identifiers.
                If False and any entry from the file has an identifier already present
                in the database, a KeyError is raised before any entries are added.
            skip_invalid: If True, entries with a non-parsable header are skipped. If
                False, a ValueError is raised when an entry is encountered which header
                could not be parsed by the header_parser. Headers of skipped entries are
                stored in the skipped_fasta_entries attribute.

        Raises:
            KeyError: If `overwrite` is `False` and a duplicate identifier is found
                across the imported files.
            ValueError: If `skip_invalid` is `False` and a FASTA header cannot be parsed
                by `header_parser`.
        """
        fasta_name = fasta_name if fasta_name is not None else pathlib.Path(path).name
        parser = get_parser(header_parser)
        parsed_protein_entries: list[DatabaseEntry] = []
        skipped_entry_headers: list[str] = []

        with open(path, "r", encoding="utf-8") as file:
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

        if not overwrite:
            duplicates = [
                e.identifier for e in parsed_protein_entries if e.identifier in self._db
            ]
            if duplicates:
                raise KeyError(
                    f"{len(duplicates)} identifier(s) from '{fasta_name}' are already "
                    f"present in the database. No entries were added."
                )

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
        if protein_entry.identifier in self._db:
            if overwrite:
                logger.warning(
                    f"Overwriting existing entry with identifier "
                    f"'{protein_entry.identifier}' in ProteinDatabase."
                )
            else:
                raise KeyError(
                    f"Identifier '{protein_entry.identifier}' already in database."
                )

        self._db[protein_entry.identifier] = protein_entry

    def filter(
        self, condition: Callable[[AbstractDatabaseEntry], bool]
    ) -> "ProteinDatabase":
        """Return a new ProteinDatabase containing only entries matching a condition.

        Note that the returned database has empty `added_fasta_files` and
        `skipped_fasta_entries` attributes, as the filtered database is not directly
        derived from any FASTA files.

        Args:
            condition: A callable that takes an `AbstractDatabaseEntry` and returns
                `True` if the entry should be included in the returned database.

        Returns:
            A new instance of `ProteinDatabase` containing only the entries for which
            `condition` returns `True`.
        """
        new_db = ProteinDatabase()
        for entry in self._db.values():
            if condition(entry):
                new_db.add_entry(entry)
        return new_db

    def write_fasta(
        self,
        path: PathLike,
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
        fasta_records: list[profasta.io.AbstractFastaRecord] = []
        if header_writer is None:
            fasta_records = list(self._db.values())
        else:
            writer = get_writer(header_writer)
            for protein_entry in self._db.values():
                header = writer.write(protein_entry)
                fasta_records.append(
                    profasta.io.FastaRecord(header, protein_entry.sequence)
                )
        file_open_mode = "a" if append else "w"
        with open(path, file_open_mode, encoding="utf-8") as file:
            profasta.io.write_fasta(file, fasta_records, line_width)

    def get(self, identifier: str, default: Any = None) -> DatabaseEntry | Any:
        """Get a protein entry by its identifier or return a default value."""
        return self._db.get(identifier, default)

    def keys(self) -> KeysView[str]:
        return self._db.keys()

    def values(self) -> ValuesView[AbstractDatabaseEntry]:
        return self._db.values()

    def items(self) -> ItemsView[str, AbstractDatabaseEntry]:
        return self._db.items()

    def __bool__(self) -> bool:
        """Return True if the database contains any entries."""
        return bool(self._db)

    def __getitem__(self, identifier) -> AbstractDatabaseEntry:
        return self._db[identifier]

    def __contains__(self, identifier) -> bool:
        return identifier in self._db

    def __iter__(self) -> Iterator[str]:
        return iter(self._db)

    def __len__(self) -> int:
        """Return the number of protein entries in the database."""
        return len(self._db)

    def __repr__(self) -> str:
        """Return a string representation of the ProteinDatabase."""
        num_entries = len(self._db)
        num_files = len(self.added_fasta_files)
        total_skipped = sum(
            len(headers) for headers in self.skipped_fasta_entries.values()
        )

        parts = [f"entries={num_entries}"]
        if num_files > 0:
            parts.append(f"files={num_files}")
        if total_skipped > 0:
            parts.append(f"skipped_entries={total_skipped}")

        return "ProteinDatabase(" + ", ".join(parts) + ")"
