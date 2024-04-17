"""FASTA file I/O.

This module provides functionality for parsing and writing FASTA files.

Classes:
    AbstractFastaRecord (Protocol): Interface for a protein record in a FASTA file.
    FastaRecord: Representation of a protein record in a FASTA file.

Functions:
    parse_fasta: Parse a FASTA file object and return a list of FastaRecords.
    write_fasta: Write a list of FastaRecords to a file object.
"""

from dataclasses import dataclass
from typing import Generator, Iterable, IO, Protocol


class AbstractFastaRecord(Protocol):
    """Abstract FASTA file record.

    Attributes:
        header: The FASTA header, not containing the starting ">" character.
        sequence: The amino acid sequence.
    """

    header: str
    sequence: str


@dataclass
class FastaRecord:
    """FASTA file record.

    Attributes:
        header: The FASTA header, not containing the starting ">" character.
        sequence: The amino acid sequence.
    """

    header: str
    sequence: str


def parse_fasta(file_object: IO[str]) -> Generator[FastaRecord, None, None]:
    """Parse a FASTA file object and yield FastaRecords.

    Lines starting with ">" are header lines, all others are sequence lines.
    Each header line is followed by one or multiple sequence lines. Multiple sequence
    lines are joined by removing new line characters into one single sequence string.

    Args:
        file_object: A file object to parse.

    Yields:
        Yields FastaRecords.
    """
    fasta_text = file_object.read()
    if not fasta_text.startswith("\n"):
        fasta_text = "\n" + fasta_text

    for block in fasta_text.split("\n>")[1:]:
        lines = block.split("\n")
        header = lines[0].strip()
        sequence = ("".join(lines[1:]).replace(" ", "").rstrip("*")).upper()
        yield FastaRecord(header, sequence)


def write_fasta(
    file_object: IO[str],
    fasta_records: Iterable[AbstractFastaRecord],
    line_width: int = -1,
):
    """Write a list of FASTA entries to a file object.

    Args:
        file_object: A file object to write to.
        fasta_records: A list of FastaRecords.
        line_width: The number of sequence characters per line, the default value is -1.
            If -1, the sequence is not split into multiple lines.
    """
    fasta_strings = [
        make_record_string(entry.header, entry.sequence, line_width)
        for entry in fasta_records
    ]
    file_object.write("\n".join(fasta_strings))
    file_object.write("\n")


def make_record_string(header: str, sequence: str, line_width: int = -1) -> str:
    """Create a FASTA string from a header and a sequence.

    Args:
        header: The header string, not containing the starting ">" character.
        sequence: The sequence string.
        line_width: The number of sequence characters per line, the default value is -1.
            If -1, the sequence is not split into multiple lines.

    Returns:
        A FASTA entry string as it appears in a FASTA file.
    """
    lines = [f">{header}"]
    if line_width == -1:
        lines.append(sequence)
    else:
        lines.extend(
            [sequence[i : i + line_width] for i in range(0, len(sequence), line_width)]
        )
    return "\n".join(lines)
