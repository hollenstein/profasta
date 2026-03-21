# Changelog

---

## Version [0.1.1]
Released: 2026-03-21

### Added
- `ProteinDatabase.write_fasta()` now automatically creates missing parent directories.
- `write_fasta()` now ensures output starts on a new line when appending to an existing file.

### Changed
- `ProteinDatabase` now inherits from `collections.abc.Mapping`.
- Default FASTA line width in `io.write_fasta()` changed to 60 characters, to ensure consistency with `db.ProteinDatabase.write_fasta()`.
- Registry operations now have standardized error semantics.
- Clarified in readme that DecoyWriter is not registered by default; refined parser docstrings.

### Fixed
- UTF-8 encoding is now enforced for all file I/O for cross-platform consistency.
- Remaining unintended public exposure of the `_db` attribute removed from `ProteinDatabase`.

---

## Version [0.1.0]
Released: 2026-03-20
 
### Added
- `ProteinDatabase` now implements dict-like interface with `__len__`, `__bool__`, and `__repr__` methods.
- `ProteinDatabase.filter()` method for condition-based entry selection.
- `ProteinDatabase.from_fasta()` convenience constructor that supports loading from multiple FASTA paths.
- `DecoyWriter` can now be instantiated with a custom decoy tag.
- `write_decoy_fasta()` function for writing decoy sequences directly to a FASTA file.
- Added a `validation` module, currently providing the function `find_header_ascii_issues` for FASTA header ASCII checks.
- Improved parser and writer registry:
  - Added `list_parsers` and `list_writers` for inspecting the content of the registy.
  - Added `replace_parser` and `replace_writer` for updating non built-in entries.
- Warning logged when overwriting existing entries in `ProteinDatabase`.
- Required and optional header fields are now documented for FASTA header parsers and writers.

### Changed
- `ProteinDatabase.db` is now a private attribute.
- Parser and writer registries are now private attributes.
- Package exports and the public API are now explicitly defined via `__all__`.

### Fixed
- `ProteinDatabase.add_fasta()` no longer leaves the database in an inconsistent state when duplicate identifiers are encountered.

### Chores
- (!) Dropped Python 3.9 / 3.10 support.
- Added 3.13 / 3.14 support.
- Switched build system from setuptools to hatchling.
- Switched CI to uv for faster installs; added linting with ruff.
- Removed nox for local testing.

---

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
