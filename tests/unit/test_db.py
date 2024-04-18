import io
import pytest

import profasta.db


class TestProteinDatabase:
    def test_add_fasta_does_not_add_entries_when_failing_to_parse_a_header(self, tmp_path):  # fmt: skip
        fasta_entries = [
            ">xx|uniprot_like_01|entry\nMKKK",
            ">not_uniprot_like_entry\nMRRR",
        ]
        fasta_path = tmp_path / "test.fasta"
        with open(fasta_path, "w") as file:
            file.write("\n".join(fasta_entries))

        protein_db = profasta.db.ProteinDatabase()
        try:
            protein_db.add_fasta(fasta_path, header_parser="uniprot_like")
        except ValueError:
            pass

        assert len(protein_db.db) == 0
        assert len(protein_db.imported_fasta_files) == 0
