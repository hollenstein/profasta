# ProFASTA
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fhollenstein%2Fprofasta%2Fmain%2Fpyproject.toml)
[![pypi](https://img.shields.io/pypi/v/profasta)](https://pypi.org/project/profasta)
[![CI](https://github.com/hollenstein/profasta/actions/workflows/ci.yml/badge.svg)](https://github.com/hollenstein/profasta/actions/workflows/ci.yml)

## Introduction

ProFASTA is a Python library for working with FASTA files containing protein records. It prioritizes simplicity while providing a practical set of features for proteomics-based mass spectrometry workflows.

Core functionality includes:

- **Parsing and writing FASTA files** via `profasta.io`
- **Structured header parsing** via a registry of built-in and user-defined parsers
- **A protein database** (`ProteinDatabase`) for managing entries loaded from one or more FASTA files
- **Decoy database generation** by sequence reversal
- **Header validation** for non-ASCII characters

ProFASTA is developed as part of the computational toolbox for the [Mass Spectrometry Facility](https://www.maxperutzlabs.ac.at/research/facilities/mass-spectrometry-facility) at the Max Perutz Labs (University of Vienna).

## Similar projects

If ProFASTA doesn't meet your requirements, consider exploring these alternative Python packages with a focus on protein-containing FASTA files:

- [fastapy](https://pypi.org/project/fastapy/) is a lightweight package with no dependencies that offers FASTA reading functionality.
- [protfasta](https://pypi.org/project/protfasta/) is another library with no dependencies that provides reading functionality along with basic validation (e.g., duplicate headers, conversion of non-canonical amino acids). The library also allows writing FASTA files with the ability to specify the sequence line length.
- [pyteomics](https://pyteomics.readthedocs.io/en/latest/index.html) is a feature-rich package that provides tools to handle various sorts of proteomics data. It provides functions for FASTA reading, automatic parsing of headers (in various formats defined at uniprot.org), writing, and generation of decoy entries. Note that pyteomics is a large package with many dependencies.

## Requirements

Python >= 3.11

ProFASTA has no dependencies beyond the Python standard library.

## Installation

Install from [PyPI](https://pypi.org/project/profasta/):

```
pip install profasta
```

## Key concepts

### FASTA parsing

The `profasta.io.parse_fasta` function reads a FASTA file and yields `FastaRecord` objects. Sequences are automatically normalized: letters are converted to uppercase, spaces are removed, and trailing `*` characters are stripped.

```python
import profasta.io

with open("proteins.fasta", "r") as f:
    for record in profasta.io.parse_fasta(f):
        print(record.header, record.sequence)
```

### Header parsers and the registry

ProFASTA uses a registry system for header parsers and writers. Built-in parsers are registered under the following names:

| Name | Description |
|---|---|
| `"default"` | Splits on the first whitespace; never fails |
| `"uniprot"` | Strict UniProt format parser |
| `"uniprot_like"` | Tolerant UniProt-like format parser |

Built-in writers follow the same naming convention and include an additional `"decoy"` writer that prepends a `rev_` tag to the header.

Custom parsers and writers can be registered via:

```python
profasta.parser.register_parser("my_parser", MyParser)
profasta.parser.register_writer("my_writer", MyWriter)
```

A parser must implement a `parse(header: str) -> ParsedHeader` classmethod, and a writer must implement a `write(parsed_header: ParsedHeader) -> str` classmethod.

### ProteinDatabase

The `ProteinDatabase` class provides a dict-like interface for managing protein entries loaded from FASTA files:

```python
import profasta

db = profasta.ProteinDatabase()
db.add_fasta("proteins.fasta", header_parser="uniprot")

entry = db["O75385"]
print(entry.header_fields["gene_name"])  # ULK1
```

Multiple FASTA files can be added to the same database. Entries with unparseable headers can be skipped using `skip_invalid=True`.

A `ProteinDatabase` can also be created directly from one or more FASTA files using the `from_fasta` convenience constructor:

```python
fasta_paths = ["proteome1.fasta", "proteome2.fasta"]
db = profasta.ProteinDatabase.from_fasta(*fasta_paths, header_parser="uniprot")
```

Entries can be filtered by a condition using the `filter` method, which returns a new `ProteinDatabase`:

```python
human_db = db.filter(lambda e: e.header_fields.get("organism_identifier") == "9606")
```

### Header validation

The `profasta.validation` module provides a function for checking FASTA records for non-ASCII characters in their headers, which can cause issues in downstream processing:

```python
import profasta.validation

with open("proteins.fasta", "r") as f:
    records = list(profasta.io.parse_fasta(f))

issues = profasta.validation.find_header_ascii_issues(records)
for issue in issues:
    print(issue.header, issue.non_ascii_characters)
```

## Usage examples

### Load a UniProt FASTA file and access a protein entry

```python
import profasta

db = profasta.ProteinDatabase()
db.add_fasta("./examples/uniprot_hsapiens_10entries.fasta", header_parser="uniprot")

entry = db["O75385"]
print(entry.header_fields["gene_name"])  # ULK1
```

### Combine multiple FASTA files and add decoy entries

A common proteomics workflow is to combine one or more FASTA files and append reversed decoy sequences. Use `profasta.write_decoy_fasta` to write decoy entries directly to a FASTA file:

```python
import profasta

# Load one or more forward databases
db = profasta.ProteinDatabase()
db.add_fasta("proteome.fasta", header_parser="uniprot")
db.add_fasta("additional.fasta", header_parser="uniprot")

# Write the forward entries, then append decoy entries with reversed sequences
output_path = "combined_with_decoys.fasta"
db.write_fasta(output_path, header_writer="default")
profasta.decoy.write_decoy_fasta(db, output_path, append=True)
```

Decoy headers are automatically prefixed with `rev_`. A custom prefix can be set via the `decoy_tag` argument:

```python
profasta.decoy.write_decoy_fasta(db, output_path, append=True, decoy_tag="decoy_")
```

## Contributors

- Juraj Ahel - [@xeniorn](https://github.com/xeniorn)