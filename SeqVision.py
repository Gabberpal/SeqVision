from typing import List, Union, Callable, Dict, Tuple
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
    fastq_data: Dict[str, Tuple[str, str]],
    gc_bounds: Union[Tuple[int, int], int, float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: Union[int, float] = 0,
) -> Dict[str, Tuple[str, str]]:
    """
    Filters reads by GC content, sequence length, and quality threshold.

    Parameters:
    ----------
    fastq_data : Dict[str, Tuple[str, str]]
        Dictionary where keys are read names,
        and values are tuples
        containing the sequence and quality.

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
        Filtered dictionary of reads
        that meet the specified conditions for GC content, length, and quality.
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    filtered_fastq_data = fastq_data.copy()
    filtered_fastq_data = filter_reads_by_gc(filtered_fastq_data, gc_bounds)
    filtered_fastq_data = filter_length_bounds(
        filtered_fastq_data, length_bounds
        )
    filtered_fastq_data = filter_quality_threshold(
        filtered_fastq_data, quality_threshold
    )
    return filtered_fastq_data


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
