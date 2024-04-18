import pytest

import profasta.decoy as decoy


class TestReverseSequence:
    @pytest.mark.parametrize(
        "sequence, expected_result",
        [
            ("ACDEF", "FEDCA"),
            ("MMMBBBCCC", "CCCBBBMMM"),
            ("A", "A"),
            ("", ""),
        ],
    )
    def test_without_keeping_nterminal_resiudes(self, sequence, expected_result):
        assert decoy._reverse_sequence(sequence, keep_nterm=False, keep_nterm_methionine=False) == expected_result  # fmt: skip

    @pytest.mark.parametrize(
        "sequence, expected_result",
        [
            ("ACDEF", "AFEDC"),
            ("MMMBBBCCC", "MCCCBBBMM"),
            ("A", "A"),
            ("", ""),
        ],
    )
    def test_with_keeping_nterminal_resiude(self, sequence, expected_result):
        assert decoy._reverse_sequence(sequence, keep_nterm=True, keep_nterm_methionine=False) == expected_result  # fmt: skip

    @pytest.mark.parametrize(
        "sequence, expected_result",
        [
            ("ACDEF", "FEDCA"),
            ("MMMBBBCCC", "MCCCBBBMM"),
            ("A", "A"),
            ("", ""),
        ],
    )
    def test_with_keeping_nterminal_methionione_resiude(self, sequence, expected_result):  # fmt: skip
        assert decoy._reverse_sequence(sequence, keep_nterm=False, keep_nterm_methionine=True) == expected_result  # fmt: skip

    def test_keep_nterminal_residue_overrules_keep_nterminal_methionine_residue(self):
        sequence = "ABCDEF"
        expected_result = "AFEDCB"
        assert decoy._reverse_sequence(sequence, keep_nterm=True, keep_nterm_methionine=True) == expected_result  # fmt: skip
