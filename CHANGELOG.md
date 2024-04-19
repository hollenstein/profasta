# Changelog

----------------------------------------------------------------------------------------

## Version [0.0.5]
Released: 2024-04-19

### Added
- Add otpion to `db.ProteinDatabase.add_fasta` that allows skipping entries which headers could not be parsed, instead of raising a `ValueError`. (Suggested by @xeniorn)
- Added `keys`, `values`, and `items` methods to `db.ProteinDatabase` to allow more convenient iteration over the database's entries.

### Changed
- Made `decoy.reverse_sequence` a private function.
- Renamed the protocol classes `HeaderParser` and `HeaderWriter` to `AbstractHeaderParser` and `AbstractHeaderWriter` to be consistent with the naming of the other abstract classes. (Suggested by @xeniorn)

### Fixed
- Parsing a FASTA file returned invalid protein sequences when the sequence contained a terminal `*` character or lowercase letters. Terminal `*` characters are now removed from the sequence and the sequence is capitalized. (Contributed by @xeniorn)

### Chores
- Added a GitHub Actions CI workflow for automated testing. (Contributed by @xeniorn)
- Minor corrections and additions to some docstrings.
- Added a Jupyter notebook containing usage examples for the ProFASTA library.
