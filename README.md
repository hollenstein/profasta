# ProFASTA
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fhollenstein%2Fprofasta%2Fmain%2Fpyproject.toml)
[![pypi](https://img.shields.io/pypi/v/profasta)](https://pypi.org/project/profasta)
[![unit-tests](https://github.com/hollenstein/profasta/actions/workflows/python-package.yml/badge.svg?branch=main)](https://github.com/hollenstein/profasta/actions/workflows/python-package.yml)

## Introduction
ProFASTA is a Python library for working with FASTA files containing protein records. Unlike other packages, ProFASTA prioritizes simplicity, while aiming to provide a set of useful features required in the field of proteomics based mass spectrometry. 

The library is still in early development and the interface might change over time. At the current stage ProFASTA provides functionality for parsing and writing FASTA files, as well as for providing access to protein records imported from FASTA files.

ProFASTA is developed as part of the computational toolbox for the [Mass Spectrometry Facility](https://www.maxperutzlabs.ac.at/research/facilities/mass-spectrometry-facility) at the Max Perutz Labs (University of Vienna).

## Similar projects
If ProFASTA doesn't meet your requirements, consider exploring these alternative Python packages with a focus on protein-containing FASTA files:

- [fastapy](https://pypi.org/project/fastapy/) is a lightweight package with no dependencies that offers FASTA reading functionality.
- [protfasta](https://pypi.org/project/protfasta/) is another library with no dependencies that provides reading functionality along with basic validation (e.g., duplicate headers, conversion of non-canonical amino acids). The library also allows writing FASTA files with the ability to specify the sequence line length.
- [pyteomics](https://pyteomics.readthedocs.io/en/latest/index.html) is a feature-rich package that provides tools to handle various sorts of proteomics data. It provides functions for FASTA reading, automatic parsing of headers (in various formats defined at uniprot.org), writing, and generation of decoy entries. Note that pyteomics is a large package with many dependencies.

## Usage example
The following code snippet shows how to import a FASTA file containing UniProt protein entries, retrieve a protein record by its UniProt accession number and print its gene name:

```python
>>> import profasta
>>> 
>>> fasta_path = "./example_data/uniprot_hsapiens_10entries.fasta"
>>> db = profasta.db.ProteinDatabase()
>>> db.add_fasta(fasta_path, header_parser="uniprot")
>>> protein_record = db["O75385"]
>>> print(protein_record.header_fields["gene_name"])
ULK1
```

## Requirements
Python >= 3.9

## Installation
The following command will install the latest version of ProFASTA and its dependencies from PyPi, the Python Packaging Index:

```
pip install profasta
```

To uninstall the ProFASTA library use:

```
pip uninstall profasta
```

## Planned features
**Main requirements**
- [x] parse FASTA file
- [x] parse FASTA header
    - [x] built-in parser that never fails
    - [x] built-in parser for uniprot format
    - [x] allow user defined parser
- [x] write FASTA file
    - [x] allow custom FASTA header generation
    
**Additional features**
- [x] read multiple FASTA files and write a combined file
- [x] add protein records to an existing FASTA file
- [x] generate decoy protein records by reversing the sequence
    - [x] add decoy protein records to an existing FASTA file
- [ ] validate FASTA file / FASTA records

## Contributors

- Juraj Ahel - [@xeniorn](https://github.com/xeniorn)