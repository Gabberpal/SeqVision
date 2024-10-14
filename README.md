# SeqVision 

**SeqVision** - is a tool for processing and analyzing nucleotide sequences such as DNA and RNA. It provides functions for filtering, transforming, and analyzing FASTQ, FASTA and BLAST sequences.

## Table of Contents 

- [Functionality](#Functionality)
- [Installation](#Installation)
- [Usage](#Usage)
- [Examples](#Examples)

## Functionality 

- **Sequence Filtering**: Filter sequences by GC content, length, and quality.
- **Sequence Transformation**: Transcription, complementarity, reverse complementarity, reverse, and palindrome selection.
- **Sequence Analysis**: Calculate statistics on sequences such as GC content and length.
- **Convert multiline fasta data to oneline**: convert fasta reads to oneline.
- **Extract and sort protein description after BLAST**: Parses a BLAST output file to extract and sort protein descriptions.

## Installation 

Clone the repository:

   ```bash
   git clone https://github.com/Gabberpal/SeqVision && cd SeqVision
   ```

## Usage 

### FASTQ sequence Filtering

To filter FASTQ sequences by GC content, length, and quality, use the `filter_fastq` function.

**Parameters**:

- input_fastq (str): Name of the input FASTQ file located in the "data" folder. 
- output_fastq (str): Name of the output FASTQ file that will be located in the "filtered" folder.
- gc_bounds (Union[Tuple[int, int], int, float], optional): GC content interval (in percent) for filtering.
- length_bounds (Union[Tuple[int, int], int], optional): Sequence length interval for filtering.
- quality_threshold (Union[int, float], optional): Average quality threshold for filtering.

Example:

```python
from SeqVision import filter_fastq

filtered_fastq = filter_fastq("example_fastq.fastq", "output_fastq.fastq", gc_bounds=(20, 80), length_bounds=(50, 75), quality_threshold=20)
```
### Sequence Transformation 

```python
sequences = ["ACGT", "UGCA"]
transformed_sequences = run_dna_rna_tools(sequences, "function_name")
print(transformed_sequences)
```

#### Functions for Sequence Transformation 

- `transcribe`: Returns transcribed sequences.
- `complement`: Returns complementary sequences.
- `reverse`: Returns reverse sequences.
- `reverse_complement`: Returns reverse-complementary sequences.
- `which_palindrome`: Returns only palindromic sequences.

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