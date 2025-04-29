# SeqVision

**SeqVision** is a powerful and extensible toolkit for processing and analyzing biological sequences (DNA, RNA, proteins) using Biopython and an object-oriented approach.

---

## 📖 Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Requirements](#requirements)
- [Usage](#usage)
  - [Core Functionality](#core-functionality)
  - [FASTQ Filtering](#fastq-filtering)
  - [File Processing](#file-processing)
- [Examples](#examples)

---

## ✨ Features

- **FASTQ Filtering**: Filter reads by GC content, sequence length, and average quality score
- **Sequence Manipulation**:
  - Transcription (DNA → RNA)
  - Complement / Reverse complement
  - Validation, slicing, and transformations
- **Protein Analysis**: Compute molecular weight of amino acid sequences
- **File Processing Utilities**:
  - Convert multiline FASTA to single-line format
  - Parse and sort protein descriptions from BLAST output
- **Extensibility**: Easily extend base classes to support new sequence types

---

## ⚙️ Installation

1. Clone the repository:
```bash
git clone https://github.com/Gabberpal/SeqVision && cd SeqVision
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

---

## 📦 Requirements

- Python ≥ 3.8
- Biopython
- `pytest` or `unittest` for testing

---

## 🚀 Usage

### Core Functionality

```python
from SeqVision import (
    DNASequence,
    RNASequence,
    AminoAcidSequence,
    filter_fastq
)

# DNA manipulation
dna = DNASequence("ATGCGTA")
print(dna.reverse_complement())  # Output: TACGCAT

# Transcription to RNA
rna = dna.transcribe()
print(rna)  # Output: AUGCGUA

# Protein molecular weight
protein = AminoAcidSequence("MAKG")
print(protein.get_molecular_weight())  # Output: 459.56
```

---

### FASTQ Filtering

```python
filter_fastq(
    input_fastq="example_fastq.fastq",
    output_fastq="filtered_output.fastq",
    gc_bounds=(35, 65),
    length_bounds=(75, 150),
    quality_threshold=25
)
```

**Parameters:**

- `input_fastq` *(str)*: Path to the input FASTQ file located in the `data/` folder
- `output_fastq` *(str)*: Name of the output file written to the `filtered/` folder
- `gc_bounds` *(tuple or number)*: GC content bounds in percent (e.g., `(30, 60)`)
- `length_bounds` *(tuple or number)*: Sequence length bounds (e.g., `(100, 150)`)
- `quality_threshold` *(int or float)*: Minimum average quality score

---

### File Processing

```python
import bio_files_processor

# Convert multiline FASTA to single-line format
bio_files_processor.convert_multiline_fasta_to_oneline(
    "data/example_multiline_fasta.fasta",
    "oneline_output.fasta"
)

# Parse BLAST output and sort descriptions
bio_files_processor.parse_blast_output(
    "data/example_blast_results.txt",
    "sorted_descriptions.txt"
)
```

---

## 🧪 Examples

### DNA Manipulation

```python
dna = DNASequence("ATGCTAGCTA")
print(dna.complement())       # Output: TACGATCGAT
print(dna[2:5])               # Output: GCT (returns new DNASequence instance)
```

### Protein Analysis

```python
protein = AminoAcidSequence("MAKGHT")
print(f"MW: {protein.get_molecular_weight():.2f} Da")  # Output: MW: 656.75 Da
```

### Quality Filtering

```python
filter_fastq(
    input_fastq="data/raw_reads.fastq",
    output_fastq="filtered/high_quality.fastq",
    gc_bounds=(40, 60),
    quality_threshold=30
)
```

### FASTA Conversion

```python
bio_files_processor.convert_multiline_fasta_to_oneline(
    "data/input_multiline.fasta",
    "filtered/oneline_output.fasta"
)
```

### BLAST Description Parsing

```python
bio_files_processor.parse_blast_output(
    "data/blast_results.txt",
    "filtered/sorted_blast_descriptions.txt"
)
```

---