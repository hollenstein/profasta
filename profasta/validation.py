"""Validation of FASTA records.

This module provides functionality for validating FASTA records. It operates
on any object that satisfies the AbstractFastaRecord protocol, making it
composable with both profasta.io.parse_fasta and profasta.db.ProteinDatabase.

Classes:
    AbstractFastaRecord (Protocol): Interface required for validation.
    HeaderAsciiIssue: Represents a non-ASCII validation issue for a single header.

Functions:
    check_ascii_issues: Check a collection of FASTA records for non-ASCII characters
        in their headers.
"""

from dataclasses import dataclass, field
from typing import Iterable, Protocol


class AbstractFastaRecord(Protocol):
    """Interface required for FASTA record validation.

    Attributes:
        header: The FASTA header, not containing the starting ">" character.
        sequence: The amino acid sequence.
    """

    header: str
    sequence: str


@dataclass
class HeaderAsciiIssue:
    """Represents a non-ASCII validation issue found in a FASTA header.

    Attributes:
        header: The full FASTA header string in which the issue was found.
        non_ascii_characters: The set of unique non-ASCII characters found in
            the header.
    """

    header: str
    non_ascii_characters: set[str] = field(default_factory=set)


def find_header_ascii_issues(
    fasta_records: Iterable[AbstractFastaRecord],
) -> list[HeaderAsciiIssue]:
    """Find FASTA headers containing non-ASCII characters.

    For each record whose header contains one or more non-ASCII characters,
    a HeaderAsciiIssue is returned containing the full header string and the
    set of unique non-ASCII characters found in it.

    Args:
        fasta_records: Any iterable of objects satisfying the AbstractFastaRecord
            protocol. Compatible with profasta.io.parse_fasta() and with iterating
            over profasta.db.ProteinDatabase.values().

    Returns:
        A list of HeaderAsciiIssue objects, one per affected header. Returns an
        empty list if no issues are found.
    """
    issues = []
    for record in fasta_records:
        # Fast path: str.isascii() is highly optimized in CPython
        if not record.header.isascii():
            non_ascii = {char for char in record.header if ord(char) > 127}
            issues.append(HeaderAsciiIssue(record.header, non_ascii))
    return issues