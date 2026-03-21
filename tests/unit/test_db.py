import pytest

import profasta.db


class TestProteinDatabase:
    def test_from_fasta_imports_entries_from_multiple_paths(self, tmp_path):
        path1 = tmp_path / "a.fasta"
        path2 = tmp_path / "b.fasta"
        path1.write_text(">xx|id_01|entry\nMKKK")
        path2.write_text(">xx|id_02|entry\nMRRR")

        db = profasta.db.ProteinDatabase.from_fasta(
            path1, path2, header_parser="uniprot_like"
        )
        assert len(db) == 2
        assert "id_01" in db and "id_02" in db
        assert "a.fasta" in db.added_fasta_files and "b.fasta" in db.added_fasta_files

    def test_filter_returns_only_matching_entries(self):
        db = profasta.db.ProteinDatabase()
        db.add_entry(profasta.db.DatabaseEntry("id_01", "header_01", "MKKK", {"gene": "BRCA1"}))  # fmt: skip
        db.add_entry(profasta.db.DatabaseEntry("id_02", "header_02", "MRRR", {"gene": "TP53"}))  # fmt: skip
        db.add_entry(profasta.db.DatabaseEntry("id_03", "header_03", "MMMM", {"gene": "BRCA1"}))  # fmt: skip

        filtered = db.filter(lambda e: e.header_fields.get("gene") == "BRCA1")
        assert len(filtered) == 2
        assert "id_01" in filtered and "id_03" in filtered
        assert "id_02" not in filtered

    def test_add_fasta_raises_value_error_when_a_header_cannot_be_parsed(self, tmp_path):  # fmt: skip
        fasta_entries = [
            ">xx|uniprot_like_01|entry\nMKKK",
            ">not_uniprot_like_entry\nMRRR",
        ]
        fasta_path = tmp_path / "test.fasta"
        with open(fasta_path, "w", encoding="utf-8") as file:
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
        with open(fasta_path, "w", encoding="utf-8") as file:
            file.write("\n".join(fasta_entries))

        protein_db = profasta.db.ProteinDatabase()
        protein_db.add_fasta(
            fasta_path, header_parser="uniprot_like", skip_invalid=True
        )
        assert len(protein_db) == 1 and "uniprot_like_01" in protein_db

    def test_add_fasta_adds_records_invalid_entry_headers_when_skip_invalid_is_true(self, tmp_path):  # fmt: skip
        fasta_entries = [
            ">xx|uniprot_like_01|entry\nMKKK",
            ">not_uniprot_like_entry\nMRRR",
        ]
        fasta_path = tmp_path / "test.fasta"
        with open(fasta_path, "w", encoding="utf-8") as file:
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
        with open(fasta_path, "w", encoding="utf-8") as file:
            file.write("\n".join(fasta_entries))

        protein_db = profasta.db.ProteinDatabase()
        try:
            protein_db.add_fasta(fasta_path, header_parser="uniprot_like")
        except ValueError:
            pass

        assert len(protein_db) == 0
        assert len(protein_db.added_fasta_files) == 0
