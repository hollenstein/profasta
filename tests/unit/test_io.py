import io

import pytest

import profasta.io


@pytest.mark.parametrize(
    "fasta_content, expected_header, expected_sequence",
    [
        (">H1\nMKKK\n>H2\nMAAA", "H1", "MKKK"),
        (">H1\nMKKK\nRRR", "H1", "MKKKRRR"),
        (">H1\nMKK K\nRR R", "H1", "MKKKRRR"),
        (">H1\nMKKK\nRRR*", "H1", "MKKKRRR"),
    ],
)
def test_parse_fasta_file(fasta_content, expected_header, expected_sequence):
    file_buffer = io.StringIO(fasta_content)
    fasta_parser = profasta.io.parse_fasta(file_buffer)
    record = next(fasta_parser)

    assert record.header == expected_header
    assert record.sequence == expected_sequence


class TestWriteFasta:
    def test_standard_write(self):
        buffer = io.StringIO()
        record_1 = profasta.io.FastaRecord(header="H1", sequence="ACGT")
        record_2 = profasta.io.FastaRecord(header="H2", sequence="TGCA")
        profasta.io.write_fasta(buffer, [record_1, record_2])

        expected = ">H1\nACGT\n>H2\nTGCA\n"
        assert buffer.getvalue() == expected

    def test_multiple_calls_to_same_buffer(self):
        buffer = io.StringIO()
        record_1 = profasta.io.FastaRecord(header="H1", sequence="ACGT")
        record_2 = profasta.io.FastaRecord(header="H2", sequence="TGCA")
        profasta.io.write_fasta(buffer, [record_1])
        profasta.io.write_fasta(buffer, [record_2])

        expected = ">H1\nACGT\n>H2\nTGCA\n"
        assert buffer.getvalue() == expected

    def test_append_to_file_adds_missing_newline(self):
        buffer = io.StringIO("EXISTING_DATA")  # No trailing newline
        buffer.seek(0, io.SEEK_END)  # Move to end as if opening in append mode
        record = profasta.io.FastaRecord(header="H1", sequence="ACGT")
        profasta.io.write_fasta(buffer, [record])
        assert buffer.getvalue() == "EXISTING_DATA\n>H1\nACGT\n"

    def test_append_to_file_with_existing_newline_does_not_add_another(self):
        buffer = io.StringIO("EXISTING_DATA\n")  # Trailing newline
        buffer.seek(0, io.SEEK_END)  # Move to end as if opening in append mode
        record = profasta.io.FastaRecord(header="H1", sequence="ACGT")
        profasta.io.write_fasta(buffer, [record])
        assert buffer.getvalue() == "EXISTING_DATA\n>H1\nACGT\n"


class TestMakeRecordString:
    @pytest.fixture
    def setup(self):
        self.header = "Header"
        self.sequence = "ACGTACGTACGTACGT"

    def test_with_no_line_breaks(self, setup):
        expected = ">Header\nACGTACGTACGTACGT"
        assert profasta.io.make_record_string(self.header, self.sequence, line_width=-1) == expected  # fmt: skip

    def test_with_line_width_5(self, setup):
        expected = ">Header\nACGTA\nCGTAC\nGTACG\nT"
        assert profasta.io.make_record_string(self.header, self.sequence, line_width=5) == expected  # fmt: skip

    def test_with_line_width_longer_than_the_sequence(self, setup):
        expected = ">Header\nACGTACGTACGTACGT"
        assert profasta.io.make_record_string(self.header, self.sequence, line_width=99) == expected  # fmt: skip
