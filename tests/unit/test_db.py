import io
import pytest

import profasta.db


class TestProteinDatabase:
    def test_add_fasta_raises_value_error_when_a_header_cannot_be_parsed(self, tmp_path):  # fmt: skip
        fasta_entries = [
            ">xx|uniprot_like_01|entry\nMKKK",
            ">not_uniprot_like_entry\nMRRR",
        ]
        fasta_path = tmp_path / "test.fasta"
        with open(fasta_path, "w") as file:
            file.write("\n".join(fasta_entries))

        protein_db = profasta.db.ProteinDatabase()
        with pytest.raises(ValueError):
            protein_db.add_fasta(fasta_path, header_parser="uniprot_like")

    def test_add_fasta_adds_valid_entries_when_skip_invalid_is_true(self, tmp_path):
        fasta_entries = [
            ">xx|uniprot_like_01|entry\nMKKK",
            ">not_uniprot_like_entry\nMRRR",
        ]
        fasta_path = tmp_path / "test.fasta"
        with open(fasta_path, "w") as file:
            file.write("\n".join(fasta_entries))

        protein_db = profasta.db.ProteinDatabase()
        protein_db.add_fasta(
            fasta_path, header_parser="uniprot_like", skip_invalid=True
        )
        assert len(protein_db.db) == 1 and "uniprot_like_01" in protein_db.db

    def test_add_fasta_adds_records_invalid_entry_headers_when_skip_invalid_is_true(self, tmp_path):  # fmt: skip
        fasta_entries = [
            ">xx|uniprot_like_01|entry\nMKKK",
            ">not_uniprot_like_entry\nMRRR",
        ]
        fasta_path = tmp_path / "test.fasta"
        with open(fasta_path, "w") as file:
            file.write("\n".join(fasta_entries))

        protein_db = profasta.db.ProteinDatabase()
        protein_db.add_fasta(fasta_path, header_parser="uniprot_like", skip_invalid=True)  # fmt: skip
        skipped_headers = protein_db.skipped_fasta_entries["test.fasta"]
        assert len(skipped_headers) == 1 and "not_uniprot_like_entry" in skipped_headers

    def test_add_fasta_does_not_add_any_entries_when_failing_to_parse_a_header(self, tmp_path):  # fmt: skip
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
        assert len(protein_db.added_fasta_files) == 0
