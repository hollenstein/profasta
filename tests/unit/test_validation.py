import pytest

import profasta.validation


class FastaRecord:
    """Minimal concrete implementation of AbstractFastaRecord for testing."""

    def __init__(self, header: str, sequence: str = ""):
        self.header = header
        self.sequence = sequence


class TestCheckAsciiIssues:
    def test_returns_empty_list_for_empty_input(self):
        assert profasta.validation.find_header_ascii_issues([]) == []

    def test_returns_empty_list_when_all_headers_are_pure_ascii(self):
        records = [
            FastaRecord("sp|P12345|GENE1_HUMAN Protein one"),
            FastaRecord("sp|Q99999|GENE2_HUMAN Protein two"),
        ]
        assert profasta.validation.find_header_ascii_issues(records) == []

    @pytest.mark.parametrize(
        "header, expected_chars",
        [
            ("Prot\u00e9in", {"\u00e9"}),
            ("Prot\u00e9in \u2013, \u00e9", {"\u00e9", "\u2013"}),
        ],
    )
    def test_detects_non_ascii_characters_in_header(self, header, expected_chars):
        records = [FastaRecord(header)]
        issues = profasta.validation.find_header_ascii_issues(records)
        assert len(issues) == 1
        assert issues[0].header == header
        assert issues[0].non_ascii_characters == expected_chars

    def test_returns_one_issue_per_affected_header(self):
        records = [
            FastaRecord("Pure ASCII header"),
            FastaRecord("Header with em\u2013dash"),
            FastaRecord("Another clean header"),
            FastaRecord("Header with \u00e9 and \u00fc"),
        ]
        issues = profasta.validation.find_header_ascii_issues(records)
        assert len(issues) == 2

    def test_issue_contains_correct_header_string(self):
        header = "sp|P12345|GENE_HUMAN Prot\u00e9in"
        records = [FastaRecord(header)]
        issues = profasta.validation.find_header_ascii_issues(records)
        assert issues[0].header == header
