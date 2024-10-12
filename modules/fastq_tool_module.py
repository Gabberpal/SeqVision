def calculate_gc_content(seq: str):
    """
    Рассчитывает содержание G и C нуклеотидов в последовательности.

    Параметры:
    ----------
    seq: str
        Строка, содержащая нуклеотидную последовательности рида.

    Возвращает:
    ----------
    int
        Процент содержания G и С нуклеотидов в последовательности.
    """
    gc_count = seq.upper().count("G") + seq.upper().count("C")
    gc_content = gc_count / len(seq) * 100
    return gc_content


def filter_reads_by_gc(fastq_data: dict, gc_bounds=(0, 100)):
    """
    Фильтрует риды по GC составу.

    Параметры:
    ----------
    fastq_dict : dict
        Словарь, где ключами являются имена ридов, а значениями — кортежи,
        содержащие последовательность и качество.

    gc_bounds : tuple or int or float, optional
        Интервал GC состава (в процентах) для фильтрации.
        Если передано одно число, оно считается верхней границей,
        а нижняя граница устанавливается в 0.
        По умолчанию равен (0, 100), т. е. все риды сохраняются.
        Примеры:
        - (20, 80) - сохраняем только риды с GC составом от 20 до 80%.
        - 44.4 - сохраняем риды с GC составом меньше, чем 44.4%.

    Возвращает:
    --------
    dict
        Отфильтрованный словарь ридов,
        удовлетворяющих заданным условиям GC состава.
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    filtered_fastq = fastq_data.copy()
    for id, (seq, quality) in fastq_data.items():
        gc_content_count = calculate_gc_content(seq)
        if not (gc_bounds[0] <= gc_content_count <= gc_bounds[1]):
            del filtered_fastq[id]
    return filtered_fastq


def calculate_length_bounds(seq: str):
    """
    Вычисляет длину рида.

    Принимает:
    ----------
    seq: int
        Строка, содержащая нуклеотидную последовательности рида.

    Возвращает:
    ----------
    int
        Длина рида.
    """
    return len(seq)


def filter_length_bounds(
    fastq_data: dict, length_bounds: tuple[tuple, int] = (0, 2**32)
):
    """
    Филтрует риды по длине.

    Принимает:
    ----------
    fastq_dict : dict
        Словарь, где ключами являются имена ридов, а значениями — кортежи,
        содержащие последовательность и качество.
    gc_bounds : tuple or int, optional
        Интервал длины рида (в количестве нуклеотидов) для фильтрации.
        Если передано одно число, оно считается верхней границей,
        а нижняя граница устанавливается в 0.
        По умолчанию равен (0, 2**32), т. е. все риды сохраняются.
        Примеры:
        - (0, 2*32) - сохраняем только риды длиной 0 до 2*32 нуклеотидов.
        - 2*10 - сохраняем риды длиной меньше, чем 2*10 нуклеотидов.
    """
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    filtered_length_bounds = fastq_data.copy()
    for id, (seq, quality) in fastq_data.items():
        length = calculate_length_bounds(seq)
        if not (length_bounds[0] <= length <= length_bounds[1]):
            del filtered_length_bounds[id]
    return filtered_length_bounds


def filter_quality_threshold(fastq_data: dict, trashhold: tuple[int, float] = 0):
    """ """
    filtered_quality_threshold = fastq_data.copy()
    for id, (seq, quality) in fastq_data.items():
        cnt = 0
        for chr in quality:
            cnt += ord(chr) - 33
        filtered_quality = cnt / len(quality)
        if not (filtered_quality >= trashhold):
            del filtered_quality_threshold[id]
    return filtered_quality_threshold
