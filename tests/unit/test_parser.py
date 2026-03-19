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


class TestWriterRegistry:
    def test_get_builtin_writer(self):
        assert profasta.parser.get_writer("default") is profasta.parser.DefaultWriter

    def test_register_and_get_custom_writer(self):
        class CustomWriter:
            @classmethod
            def write(cls, parsed_header):
                return parsed_header.header

        profasta.parser.register_writer("custom_writer", CustomWriter)
        assert profasta.parser.get_writer("custom_writer") is CustomWriter
