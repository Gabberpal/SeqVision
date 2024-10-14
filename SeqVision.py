from typing import List, Union, Callable, Dict, Tuple
import os
from modules.fastq_tool_module import (
    calculate_gc_content,
    filter_reads_by_gc,
    filter_length_bounds,
    filter_quality_threshold,
)
from modules.dna_rna_tools_module import (
    is_rna,
    is_dna,
    complement,
    reverse,
    transcribe,
    reverse_complement,
    which_palindrome,
)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[int, int], int, float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: Union[int, float] = 0,
) -> Dict[str, Tuple[str, str]]:
    """
    Filters FASTQ reads by GC content, sequence length, and quality threshold.

    Parameters:
    ----------
    input_fastq : str
        Name of the FASTQ file located in the "data: folder.

    output_fastq : str
        Name of the output FASTQ file 
        that will be located in the "filtered" folder.

    gc_bounds : Union[Tuple[int, int], int, float], optional
        GC content interval (in percent) for filtering.
        If a single number is provided, it is considered the upper bound,
        and the lower bound is set to 0.
        Defaults to (0, 100), i.e., all reads are preserved.

    length_bounds : Union[Tuple[int, int], int], optional
        Sequence length interval for filtering.
        If a single number is provided, it is considered the upper bound,
        and the lower bound is set to 0.
        Defaults to (0, 2**32), i.e., all reads are preserved.

    quality_threshold : Union[int, float], optional
        Average quality threshold for filtering.
        Defaults to 0, i.e., all reads are preserved.

    Returns:
    --------
    Dict[str, Tuple[str, str]]
        Filtered dictionary of reads that meet the specified conditions
        for GC content, length, and quality.
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    file_path = os.path.join("data", input_fastq)
    output_file_path = os.path.join("filtered", output_fastq)
    with open(file_path, "r") as fastq_file:
        with open(output_file_path, "a") as output_file:
            while True:
                id_line = fastq_file.readline().strip()
                if not id_line:
                    break
                seq_line = fastq_file.readline().strip()
                id2_line = fastq_file.readline().strip()
                quality_line = fastq_file.readline().strip()
                if not (id_line.startswith("@SR") 
                        and id2_line.startswith("+SR")):
                    raise TypeError("The fuction works only with .fastq files")
                if (
                    filter_reads_by_gc(seq_line, gc_bounds) and
                    filter_length_bounds(seq_line, length_bounds) and
                    filter_quality_threshold(quality_line, quality_threshold)
                ):
                    output_file.write(id_line + "\n")
                    output_file.write(seq_line + "\n")
                    output_file.write(id2_line + "\n")
                    output_file.write(quality_line + "\n")


def run_dna_rna_tools(
        *args: Union[str, Callable[[str], str]]
        ) -> Union[str, List[str]]:
    """
    Performs the specified function on DNA or RNA sequences.

    Parameters:
    ----------
    *args : Union[str, Callable[[str], str]]
        Variable number of arguments,
        where the last argument is the function name,
        and all preceding arguments are DNA or RNA sequences.

    Returns:
    --------
    Union[str, List[str]]
        The result of the function applied to the sequences.
        If a single sequence is provided, a string is returned.
        If multiple sequences are provided, a list of strings is returned.

    Exceptions:
    -----------
    ValueError
        If fewer than two arguments are provided or if the last
        argument is not a valid function name.
    ValueError
        If the sequences are not DNA or RNA.
    """
    if len(args) < 2:
        raise ValueError("The function requires at least 2 arguments.")
    function = args[-1]
    seqs = args[:-1]
    if function not in [
        "transcribe",
        "complement",
        "reverse_complement",
        "reverse",
        "which_palindrome",
    ]:
        raise ValueError("Such a function does not exist, please choose another one.")
    returned_seqs = []
    for seq in seqs:
        if not is_dna(seq) and not is_rna(seq):
            raise ValueError(
                "The function only works with DNA and RNA sequences."
                )
        returned_seqss = functions[function](seq)
        returned_seqs.append(returned_seqss)
    returned_seqs = [i for i in returned_seqs if i is None]
    if len(returned_seqs) > 1:
        return returned_seqs
    else:
        return returned_seqs[0]
    