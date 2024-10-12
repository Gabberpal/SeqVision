def calculate_gc_content(seq: str):
    """
    Calculates the content of G and C nucleotides in the sequence.

    Parameters:
    ----------
    seq: str
        String containing the nucleotide sequence of the read.

    Returns:
    ----------
    int
        Percentage of G and C nucleotides in the sequence.
    """
    gc_count = seq.upper().count("G") + seq.upper().count("C")
    gc_content = gc_count / len(seq) * 100
    return gc_content


def filter_reads_by_gc(fastq_data: dict, gc_bounds=(0, 100)):
    """
    Filters reads by GC content.

    Parameters:
    ----------
    fastq_dict : dict
        Dictionary where keys are read names and values are tuples
        containing the sequence and quality.

    gc_bounds : tuple or int or float, optional
        GC content interval (in percent) for filtering.
        If a single number is provided, it is considered the upper bound,
        and the lower bound is set to 0.
        Defaults to (0, 100), i.e., all reads are preserved.
        Examples:
        - (20, 80) - preserves only reads with GC content between 20% and 80%.
        - 44.4 - preserves reads with GC content less than 44.4%.

    Returns:
    --------
    dict
        Filtered dictionary of reads
        that meet the specified GC content conditions.
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    filtered_fastq = fastq_data.copy()
    for id, (seq, quality) in fastq_data.items():
        gc_content_count = calculate_gc_content(seq)
        if not (gc_bounds[0] <= gc_content_count <= gc_bounds[1]):
            del filtered_fastq[id]
    return filtered_fastq


def filter_length_bounds(
    fastq_data: dict, length_bounds: tuple[tuple, int] = (0, 2**32)
):
    """
    Filters reads by length.

    Parameters:
    ----------
    fastq_dict : dict
        Dictionary where keys are read names and values are tuples
        containing the sequence and quality.

    length_bo   unds : tuple or int, optional
        Read length interval (in number of nucleotides) for filtering.
        If a single number is provided, it is considered the upper bound,
        and the lower bound is set to 0.
        Defaults to (0, 2**32), i.e., all reads are preserved.
        Examples:
        - (0, 2*32) - preserves only reads with lengths between 0 and 2*32 nucleotides.
        - 2*10 - preserves reads with lengths less than 2*10 nucleotides.
    """
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    filtered_length_bounds = fastq_data.copy()
    for id, (seq, quality) in fastq_data.items():
        length = length(seq)
        if not (length_bounds[0] <= length <= length_bounds[1]):
            del filtered_length_bounds[id]
    return filtered_length_bounds


def filter_quality_threshold(fastq_data: dict, trashhold: tuple[int, float] = 0):
    """ """
    filtered_quality_threshold = fastq_data.copy()
    for id, (seq, quality) in fastq_data.items():
        filtered_quality = 0
        filtered_quality = sum([ord(chr) - 33 for chr in quality]) / len(quality)
        if not (filtered_quality >= trashhold):
            del filtered_quality_threshold[id]
    return filtered_quality_threshold
