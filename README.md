# SeqVision 

**SeqVision** - is a tool for processing and analyzing nucleotide sequences such as DNA and RNA. It provides functions for filtering, transforming, and analyzing sequences.

## Table of Contents 

- [Functionality](#Functionality)
- [Installation](#Installation)
- [Usage](#Usage)
- [Examples](#Examples)

## Functionality 

- **Sequence Filtering**: Filter sequences by GC content, length, and quality.
- **Sequence Transformation**: Transcription, complementarity, reverse complementarity, reverse, and palindrome selection.
- **Sequence Analysis**: Calculate statistics on sequences such as GC content and length.

## Installation 

Clone the repository:

   ```bash
   git clone https://github.com/Gabberpal/SeqVision && cd SeqVision
   ```

## Usage 

### Sequence Filtering

To filter sequences by GC content, length, and quality, use the `filter_fastq` function. Example:

```python
from seqvision import filter_fastq

fastq_data = {
    '@SRX079801': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079802': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
    '@SRX079803': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD')
}

filtered_fastq = filter_fastq(fastq_data, gc_bounds=(20, 80), length_bounds=(50, 75), quality_threshold=20)
print(filtered_fastq)
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


