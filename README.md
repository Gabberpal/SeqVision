# SeqVision 

**SeqVision** - A comprehensive toolkit for processing and analyzing biological sequences (DNA, RNA, proteins) with Biopython integration and object-oriented design.

## Table of Contents 
- [Features](#features)
- [Installation](#installation)
- [Requirements](#requirements)
- [Usage](#usage)
- [Examples](#examples)
- [Data Structure](#data-structure)

## Features

- **FASTQ Filtering**: Advanced filtering by GC content, sequence length, and quality scores using Biopython
- **Sequence Manipulation**: OOP-based operations for nucleic acids:
  - Transcription (DNA → RNA)
  - Complement/reverse complement
  - Sequence validation and transformation
- **Protein Analysis**: Molecular weight calculation for amino acid sequences
- **File Processing**:
  - Multiline FASTA to single-line conversion
  - BLAST results parsing and sorting
- **Extensible Core**: Base classes for custom sequence types development

## Installation

1. Clone the repository:
```bash
git clone https://github.com/Gabberpal/SeqVision && cd SeqVision
```
2. Install requirements:
```bash
pip install -r requirements.txt
```
## Usage 

### Core Functionality

```python
from SeqVision import (
    DNASequence,
    RNASequence,
    AminoAcidSequence,
    filter_fastq
)

# Create validated DNA sequence
dna = DNASequence("ATGCGTA")
print(dna.reverse_complement())  # TACGCAT

# Transcribe to RNA
rna = dna.transcribe()
print(rna)  # AUGCGUA

# Protein analysis
protein = AminoAcidSequence("MAKG")
print(protein.get_molecular_weight())  # 459.56
```

### FASTQ sequence Filtering

To filter FASTQ sequences by GC content, length, and quality, use the `filter_fastq` function.

```python
filter_fastq(
    input_fastq="input.fastq",
    output_fastq="filtered.fastq",
    gc_bounds=(35, 65),
    length_bounds=(75, 150),
    quality_threshold=25
)
```
**Parameters**:

- input_fastq (str): Name of the input FASTQ file located in the "data" folder. 
- output_fastq (str): Name of the output FASTQ file that will be located in the "filtered" folder.
- gc_bounds (Union[Tuple[int, int], int, float], optional): GC content interval (in percent) for filtering.
- length_bounds (Union[Tuple[int, int], int], optional): Sequence length interval for filtering.
- quality_threshold (Union[int, float], optional): Average quality threshold for filtering.

### File Processing

```python
import bio_files_processor

# Convert FASTA format
bio_files_processor.convert_multiline_fasta_to_oneline(
    "input.fasta", 
    "singleline.fasta"
)

# Process BLAST results
bio_files_processor.parse_blast_output(
    "blast_results.txt",
    "sorted_descriptions.txt"
)
```

## Examples

### DNA Operations

```python
dna = DNASequence("ATGCTAGCTA")
print(dna.complement())  # TACGATCGAT
print(dna[2:5])  # GCT (returns new DNASequence instance)
```

### Protein Analysis

```python
protein = AminoAcidSequence("MAKGHT")
print(f"Molecular weight: {protein.get_molecular_weight():.2f} Da")
```

### Quality Filtering

```python
# Filter with 40-60% GC content and Q30+ quality
filter_fastq(
    "raw_reads.fastq",
    "high_quality.fastq",
    gc_bounds=(40, 60),
    quality_threshold=30
)
```

### Convert multiline FASTA data to oneline

To convert multiline FASTA data to oneline use `convert_multiline_fasta_to_oneline()` function.
Ecample:

```Python
import bio_files_processor 

convert_multiline_fasta_to_oneline('example_miltiline_fasta.fasta', 'output_fasta.fasta')
```

### Extract and sort protein description after BLAST

To extract and sort protein description after BLAST use `parse_blast_output()` function.
Example:

```Python
import bio_files_processor 

parse_blast_output("example_blast_results.txt", "blast_output.txt")
```

## Data structure 

SeqVision/
├── data/                   # Input files
│   ├── example.fastq
│   ├── multiline.fasta
│   └── blast_results.txt
├── filtered/               # Processed outputs
├── SeqVision.py            # Main functionality
├── bio_files_processor.py  # File format handlers
├── requirements.txt
└── README.md

