from pathlib import Path

import pytest

import profasta

FASTA_PATH = Path(__file__).parent / "uniprot_hsapiens_10entries.fasta"


@pytest.fixture(autouse=True)
def reset_registries():
    profasta.parser._reset_registries()


class CustomParser:
    @classmethod
    def parse(cls, header):
        return profasta.parser.ParsedHeader(header, header)


class CustomWriter:
    @classmethod
    def write(cls, parsed_header):
        return parsed_header.header


def test_built_in_uniprot_header_parser_read_write_header_roundtrip():
    UniprotParser = profasta.parser.get_parser("uniprot")
    UniprotWriter = profasta.parser.get_writer("uniprot")

    with open(FASTA_PATH, "r", encoding="utf-8") as file:
        for fasta_record in profasta.io.parse_fasta(file):
            parsed_header = UniprotParser.parse(fasta_record.header)
            written_header = UniprotWriter.write(parsed_header)
            assert fasta_record.header == written_header


def test_custom_header_parser_read_write_header_roundtrip():
    profasta.parser.register_parser("test_parser", CustomParser)
    profasta.parser.register_writer("test_writer", CustomWriter)

    TestParser = profasta.parser.get_parser("test_parser")
    TestWriter = profasta.parser.get_writer("test_writer")

    with open(FASTA_PATH, "r", encoding="utf-8") as file:
        for fasta_record in profasta.io.parse_fasta(file):
            parsed_header = TestParser.parse(fasta_record.header)
            written_header = TestWriter.write(parsed_header)
            assert fasta_record.header == written_header


def test_protein_database_read_write_fasta_roundtrip(tmp_path):
    temp_fasta_path = tmp_path / "written_from_db.fasta"
    db = profasta.ProteinDatabase()
    db.add_fasta(FASTA_PATH, header_parser="default")
    db.write_fasta(temp_fasta_path, header_writer="default")

    db2 = profasta.ProteinDatabase()
    db2.add_fasta(temp_fasta_path, header_parser="default")

    for identifier in db:
        assert db[identifier].header == db2[identifier].header
        assert db[identifier].sequence == db2[identifier].sequence


def test_protein_database_write_fasta_in_append_mode(tmp_path):
    temp_fasta_path = tmp_path / "written_from_db.fasta"
    db_write = profasta.ProteinDatabase()
    db_write.add_fasta(FASTA_PATH, header_parser="default")
    db_write.write_fasta(temp_fasta_path)

    db_append = profasta.ProteinDatabase()
    db_append.add_entry(profasta.DatabaseEntry("test_header", "test_header", "seq", {}))  # fmt: skip
    db_append.write_fasta(temp_fasta_path, append=True)

    db_read = profasta.ProteinDatabase()
    db_read.add_fasta(temp_fasta_path, header_parser="default")

    for identifier in [*db_write, *db_append]:
        assert identifier in db_read


def test_create_decoy_database(tmp_path):
    decoy_fasta_path = tmp_path / "decoy.fasta"
    db = profasta.ProteinDatabase()
    db.add_fasta(FASTA_PATH, header_parser="uniprot")
    decoy_db = profasta.decoy.create_decoy_db(db, keep_nterm=False, keep_nterm_methionine=False)  # fmt: skip
    decoy_db.write_fasta(decoy_fasta_path)

    db2 = profasta.ProteinDatabase()
    db2.add_fasta(decoy_fasta_path, header_parser="uniprot")

    for identifier in db:
        assert db[identifier].header == db2[identifier].header
        assert db[identifier].sequence == db2[identifier].sequence[::-1]


def test_write_decoy_fasta(tmp_path):
    decoy_fasta_path = tmp_path / "subdir" / "decoy.fasta"
    decoy_tag = "DECOY_"

    db = profasta.ProteinDatabase()
    db.add_fasta(FASTA_PATH, header_parser="default")
    profasta.decoy.write_decoy_fasta(
        db,
        decoy_fasta_path,
        decoy_tag=decoy_tag,
        keep_nterm=False,
        keep_nterm_methionine=False,
    )

    decoy_records = {}
    with open(decoy_fasta_path, "r", encoding="utf-8") as file:
        for record in profasta.io.parse_fasta(file):
            decoy_records[record.header] = record

    assert len(decoy_records) == len(db)
    for entry in db.values():
        expected_header = decoy_tag + entry.header
        assert expected_header in decoy_records
        assert decoy_records[expected_header].sequence == entry.sequence[::-1]
