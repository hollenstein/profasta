"""Functions for creating decoy databases.

Functions:
    create_decoy_db: Create a decoy database by reversing the sequences of the input
        database records.
"""

from copy import deepcopy
from profasta.db import ProteinDatabase


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
