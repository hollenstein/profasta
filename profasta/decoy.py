"""Functions for creating decoy databases.

Functions:
    create_decoy_db: Create a decoy database by reversing the sequences of the input
        database records.
    write_decoy_fasta: Write decoy sequences derived from a ProteinDatabase directly to
        a FASTA file.
"""

import pathlib
from copy import deepcopy
from os import PathLike

import profasta.io
from profasta.db import ProteinDatabase
from profasta.parser import DecoyWriter


def create_decoy_db(
    db: ProteinDatabase, keep_nterm: bool = False, keep_nterm_methionine: bool = True
) -> ProteinDatabase:
    """Create a decoy database by reversing the sequences of the input database records.

    Args:
        db: The input database.
        keep_nterm: If True, keep the N-terminal residue in the same position. Default
            is False.
        keep_nterm_methionine: If True, keep the N-terminal residue in the same position
            if it is a methionine. Default is True.

    Returns:
        The decoy database.
    """
    decoy_db = ProteinDatabase()
    for protein in db:
        decoy_entry = deepcopy(db[protein])
        decoy_entry.sequence = _reverse_sequence(
            decoy_entry.sequence,
            keep_nterm=keep_nterm,
            keep_nterm_methionine=keep_nterm_methionine,
        )
        decoy_db.add_entry(decoy_entry)
    return decoy_db


def write_decoy_fasta(
    db: ProteinDatabase,
    path: PathLike,
    append: bool = False,
    decoy_tag: str = "rev_",
    keep_nterm: bool = False,
    keep_nterm_methionine: bool = True,
    line_width: int = 60,
):
    """Write decoy sequences derived from a ProteinDatabase to a FASTA file.

    For each entry in the database the sequence is reversed (optionally keeping the
    N-terminal residue in place) and the header is prefixed with a decoy tag. The
    resulting records are written directly to a FASTA file without creating an
    intermediate ProteinDatabase.

    Parent directories of `path` are created automatically if they do not exist.

    Args:
        db: The protein database whose entries are used as targets for decoy generation.
        path: The path to the output FASTA file.
        append: If False, the file is created or overwritten. If True, the decoy entries
            are appended to an existing file. The default value is False.
        decoy_tag: The string prepended to each FASTA header to mark it as a decoy. The
            default value is "rev_".
        keep_nterm: If True, keep the N-terminal residue in the same position during
            sequence reversal. Default is False.
        keep_nterm_methionine: If True, keep the N-terminal methionine in the same
            position during sequence reversal. Default is True.
        line_width: The number of sequence characters per line. The default value is 60.
            If -1, the sequence is not split into multiple lines.
    """
    temp_writer = DecoyWriter.with_tag(decoy_tag)
    fasta_records = []
    for entry in db.values():
        reversed_sequence = _reverse_sequence(
            entry.sequence,
            keep_nterm=keep_nterm,
            keep_nterm_methionine=keep_nterm_methionine,
        )
        header = temp_writer.write(entry)
        fasta_records.append(profasta.io.FastaRecord(header, reversed_sequence))

    output_path = pathlib.Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    file_open_mode = "a" if append else "w"
    with open(output_path, file_open_mode) as file:
        profasta.io.write_fasta(file, fasta_records, line_width)


def _reverse_sequence(
    sequence: str, keep_nterm: bool = False, keep_nterm_methionine: bool = True
) -> str:
    """Create a decoy sequence by reversing the input sequence.

    Args:
        sequence: The input sequence.
        keep_nterm: If True, keep the N-terminal residue in the same position. Default
            is False.
        keep_nterm_methionine: If True, keep the N-terminal residue in the same position
            if it is a methionine. Default is True.
    """
    if not sequence:
        return sequence
    if keep_nterm or (keep_nterm_methionine and sequence[0] == "M"):
        return sequence[0] + sequence[1:][::-1]
    return sequence[::-1]
