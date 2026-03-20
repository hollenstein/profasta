import pytest

import profasta.parser


class TestDefaultParser:
    def test_parse_header_with_description(self):
        header = "P12345 Some protein description"
        parsed = profasta.parser.DefaultParser.parse(header)
        assert parsed.identifier == "P12345"
        assert parsed.header == header
        assert parsed.header_fields == {"description": "Some protein description"}

    def test_parse_header_without_description(self):
        header = "P12345"
        parsed = profasta.parser.DefaultParser.parse(header)
        assert parsed.identifier == "P12345"
        assert parsed.header_fields == {}


class TestDefaultWriter:
    def test_write_returns_original_header(self):
        header = "P12345 Some protein description"
        parsed = profasta.parser.ParsedHeader("P12345", header, {})
        assert profasta.parser.DefaultWriter.write(parsed) == header

    def test_parse_write_roundtrip(self):
        header = "ACC_001 A multi-word description"
        parsed = profasta.parser.DefaultParser.parse(header)
        assert profasta.parser.DefaultWriter.write(parsed) == header


class TestUniprotLikeParser:
    def test_parse_full_uniprot_header(self):
        header = "sp|O75385|ULK1_HUMAN Serine/threonine-protein kinase ULK1 OS=Homo sapiens OX=9606 GN=ULK1 PE=1 SV=2"
        expected_fields = {
            "db": "sp",
            "identifier": "O75385",
            "entry_name": "ULK1_HUMAN",
            "protein_name": "Serine/threonine-protein kinase ULK1",
            "organism_name": "Homo sapiens",
            "organism_identifier": "9606",
            "gene_name": "ULK1",
            "protein_existence": "1",
            "sequence_version": "2",
        }
        parsed_header = profasta.parser.UniprotLikeParser.parse(header)
        assert parsed_header.header_fields == expected_fields

    def test_parse_partial_uniprot_header_with_some_tag_fields_missing(self):
        header = (
            "sp|O75385|ULK1_HUMAN Serine/threonine-protein kinase ULK1 OX=9606 GN=ULK1"
        )
        expected_fields = {
            "db": "sp",
            "identifier": "O75385",
            "entry_name": "ULK1_HUMAN",
            "protein_name": "Serine/threonine-protein kinase ULK1",
            "organism_identifier": "9606",
            "gene_name": "ULK1",
        }
        parsed_header = profasta.parser.UniprotLikeParser.parse(header)
        assert parsed_header.header_fields == expected_fields

    def test_parse_partial_uniprot_header_with_protein_name_missing(self):
        header = "sp|O75385|ULK1_HUMAN OS=Homo sapiens OX=9606 GN=ULK1 PE=1 SV=2"
        expected_fields = {
            "db": "sp",
            "identifier": "O75385",
            "entry_name": "ULK1_HUMAN",
            "organism_name": "Homo sapiens",
            "organism_identifier": "9606",
            "gene_name": "ULK1",
            "protein_existence": "1",
            "sequence_version": "2",
        }
        parsed_header = profasta.parser.UniprotLikeParser.parse(header)
        assert parsed_header.header_fields == expected_fields

    def test_parse_partial_uniprot_header_with_no_tag_fields(self):
        header = "sp|O75385|ULK1_HUMAN Serine/threonine-protein kinase ULK1"
        expected_fields = {
            "db": "sp",
            "identifier": "O75385",
            "entry_name": "ULK1_HUMAN",
            "protein_name": "Serine/threonine-protein kinase ULK1",
        }

        parsed_header = profasta.parser.UniprotLikeParser.parse(header)
        assert parsed_header.header_fields == expected_fields

    def test_parse_minimal_uniprot_header_with_no_description(self):
        header = "sp|O75385|ULK1_HUMAN"
        expected_fields = {
            "db": "sp",
            "identifier": "O75385",
            "entry_name": "ULK1_HUMAN",
        }

        parsed_header = profasta.parser.UniprotLikeParser.parse(header)
        assert parsed_header.header_fields == expected_fields


class TestDecoyWriter:
    def test_decoy_tag_added(self):
        header = "sp|O75385|ULK1_HUMAN"
        parsed_header = profasta.parser.ParsedHeader("", header, {})
        decoy_header = profasta.parser.DecoyWriter.write(parsed_header)
        assert decoy_header == "rev_" + header

    def test_with_custom_tag(self):
        header = "sp|O75385|ULK1_HUMAN"
        parsed_header = profasta.parser.ParsedHeader("", header, {})
        CustomDecoyWriter = profasta.parser.DecoyWriter.with_tag("decoy_")
        assert CustomDecoyWriter.write(parsed_header) == "decoy_" + header

    def test_init_with_custom_tag_does_not_affect_default_class(self):
        header = "sp|O75385|ULK1_HUMAN"
        parsed_header = profasta.parser.ParsedHeader("", header, {})
        profasta.parser.DecoyWriter.with_tag("decoy_")
        assert profasta.parser.DecoyWriter.write(parsed_header) == "rev_" + header


class TestParserRegistry:
    def test_get_builtin_parser(self):
        assert profasta.parser.get_parser("default") is profasta.parser.DefaultParser

    def test_get_unknown_parser_raises_key_error(self):
        with pytest.raises(KeyError, match="unknown_parser"):
            profasta.parser.get_parser("unknown_parser")

    def test_register_and_get_custom_parser(self):
        class MyParser:
            @classmethod
            def parse(cls, header):
                return profasta.parser.ParsedHeader(header, header)

        profasta.parser.register_parser("my_parser", MyParser)
        assert profasta.parser.get_parser("my_parser") is MyParser

    def test_register_builtin_name_raises_value_error(self):
        with pytest.raises(ValueError, match="built-in"):
            profasta.parser.register_parser("default", profasta.parser.DefaultParser)

    def test_register_duplicate_name_raises_value_error(self):
        class AnotherParser:
            @classmethod
            def parse(cls, header):
                return profasta.parser.ParsedHeader(header, header)

        profasta.parser.register_parser("dup_parser", AnotherParser)
        with pytest.raises(ValueError, match="dup_parser"):
            profasta.parser.register_parser("dup_parser", AnotherParser)

    def test_replace_parser_updates_registry(self):
        class ParserV1:
            @classmethod
            def parse(cls, header):
                return profasta.parser.ParsedHeader(header, header)

        class ParserV2:
            @classmethod
            def parse(cls, header):
                return profasta.parser.ParsedHeader(header, header)

        profasta.parser.register_parser("replaceable_parser", ParserV1)
        profasta.parser.replace_parser("replaceable_parser", ParserV2)
        assert profasta.parser.get_parser("replaceable_parser") is ParserV2

    def test_replace_builtin_parser_raises_key_error(self):
        with pytest.raises(KeyError, match="built-in"):
            profasta.parser.replace_parser("uniprot", profasta.parser.DefaultParser)

    def test_replace_unregistered_parser_raises_key_error(self):
        with pytest.raises(KeyError, match="nonexistent_parser"):
            profasta.parser.replace_parser(
                "nonexistent_parser", profasta.parser.DefaultParser
            )

    def test_list_parsers_contains_builtins(self):
        names = profasta.parser.list_parsers()
        assert "default" in names
        assert "uniprot" in names
        assert "uniprot_like" in names

    def test_list_parsers_returns_list(self):
        assert isinstance(profasta.parser.list_parsers(), list)


class TestWriterRegistry:
    def test_get_builtin_writer(self):
        assert profasta.parser.get_writer("default") is profasta.parser.DefaultWriter

    def test_get_unknown_writer_raises_key_error(self):
        with pytest.raises(KeyError, match="unknown_writer"):
            profasta.parser.get_writer("unknown_writer")

    def test_register_and_get_custom_writer(self):
        class CustomWriter:
            @classmethod
            def write(cls, parsed_header):
                return parsed_header.header

        profasta.parser.register_writer("custom_writer", CustomWriter)
        assert profasta.parser.get_writer("custom_writer") is CustomWriter

    def test_register_builtin_name_raises_value_error(self):
        with pytest.raises(ValueError, match="built-in"):
            profasta.parser.register_writer("default", profasta.parser.DefaultWriter)

    def test_register_duplicate_name_raises_value_error(self):
        class AnotherWriter:
            @classmethod
            def write(cls, parsed_header):
                return parsed_header.header

        profasta.parser.register_writer("dup_writer", AnotherWriter)
        with pytest.raises(ValueError, match="dup_writer"):
            profasta.parser.register_writer("dup_writer", AnotherWriter)

    def test_replace_writer_updates_registry(self):
        class WriterV1:
            @classmethod
            def write(cls, parsed_header):
                return parsed_header.header

        class WriterV2:
            @classmethod
            def write(cls, parsed_header):
                return parsed_header.header

        profasta.parser.register_writer("replaceable_writer", WriterV1)
        profasta.parser.replace_writer("replaceable_writer", WriterV2)
        assert profasta.parser.get_writer("replaceable_writer") is WriterV2

    def test_replace_builtin_writer_raises_key_error(self):
        with pytest.raises(KeyError, match="built-in"):
            profasta.parser.replace_writer("default", profasta.parser.DefaultWriter)

    def test_replace_unregistered_writer_raises_key_error(self):
        with pytest.raises(KeyError, match="nonexistent_writer"):
            profasta.parser.replace_writer(
                "nonexistent_writer", profasta.parser.DefaultWriter
            )

    def test_list_writers_contains_builtins(self):
        names = profasta.parser.list_writers()
        assert "default" in names
        assert "uniprot" in names
        assert "uniprot_like" in names

    def test_list_writers_returns_list(self):
        assert isinstance(profasta.parser.list_writers(), list)
