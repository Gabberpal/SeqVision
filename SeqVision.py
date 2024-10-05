from typing import List, Union, Callable, Dict, Tuple
import typing
from modules.fastq_tool_module import (
    calculate_gc_content,
    filter_reads_by_gc,
    calculate_length_bounds,
    filter_length_bounds,
    filter_quality_threshold,
)
from modules.dna_rna_tools_module import (
    set_rna,
    set_dna,
    ComplementDictDna,
    ComplementDictRna,
    TranscribeDictDna,
    TranscribeDictRna,
    is_rna,
    is_dna,
    complement,
    reverse,
    transcribe,
    reverse_complement,
    which_palindrome,
    functions,
)


def filter_fastq(
    fastq_data: Dict[str, Tuple[str, str]],
    gc_bounds: Union[Tuple[int, int], int, float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: Union[int, float] = 0,
) -> Dict[str, Tuple[str, str]]:
    """
    Фильтрует риды по GC составу, длине последовательности и порогу качества.

    Параметры:
    ----------
    fastq_data : Dict[str, Tuple[str, str]]
        Словарь, где ключами являются имена ридов, а значениями — кортежи,
        содержащие последовательность и качество.
    gc_bounds : Union[Tuple[int, int], int, float], optional
        Интервал GC состава (в процентах) для фильтрации.
        Если передано одно число, оно считается верхней границей, а нижняя граница устанавливается в 0.
        По умолчанию равен (0, 100), т. е. все риды сохраняются.
    length_bounds : Union[Tuple[int, int], int], optional
        Интервал длины последовательности для фильтрации.
        Если передано одно число, оно считается верхней границей, а нижняя граница устанавливается в 0.
        По умолчанию равен (0, 2**32), т. е. все риды сохраняются.
    quality_threshold : Union[int, float], optional
        Порог среднего качества для фильтрации.
        По умолчанию равен 0, т. е. все риды сохраняются.

    Возвращает:
    --------
    Dict[str, Tuple[str, str]]
        Отфильтрованный словарь ридов, удовлетворяющих заданным условиям GC состава, длины и качества.
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    filtered_fastq_data = fastq_data.copy()
    filtered_fastq_data = filter_reads_by_gc(filtered_fastq_data, gc_bounds)
    filtered_fastq_data = filter_length_bounds(filtered_fastq_data, length_bounds)
    filtered_fastq_data = filter_quality_threshold(
        filtered_fastq_data, quality_threshold
    )
    return filtered_fastq_data


def run_dna_rna_tools(*args: Union[str, Callable[[str], str]]) -> Union[str, List[str]]:
    """
    Выполняет указанную функцию над последовательностями ДНК или РНК.

    Параметры:
    ----------
    *args : Union[str, Callable[[str], str]]
        Переменное количество аргументов, где последний аргумент - это имя функции,
        а все предыдущие аргументы - это последовательности ДНК или РНК.

    Возвращает:
    --------
    Union[str, List[str]]
        Результат выполнения функции над последовательностями. Если передана одна последовательность,
        возвращается строка. Если передано несколько последовательностей, возвращается список строк.

    Исключения:
    -----------
    ValueError
        Если передано менее двух аргументов или если последний аргумент не является допустимым именем функции.
    ValueError
        Если последовательности не являются ДНК или РНК.
    """
    if len(args) < 2:
        raise ValueError("Функция принимает минимум 2 аргумента")
    function = args[-1]
    seqs = args[:-1]
    if function not in [
        "transcribe",
        "complement",
        "reverse_complement",
        "reverse",
        "which_palindrome",
    ]:
        raise ValueError("Такой функции нет, выберите другую")
    returned_seqs = []
    for seq in seqs:
        if not is_dna(seq) and not is_rna(seq):
            raise ValueError("Функция работает только с ДНК и РНК последовательностями")
        returned_seqss = functions[function](seq)
        returned_seqs.append(returned_seqss)
    returned_seqs = [i for i in returned_seqs if i != None]
    if len(returned_seqs) > 1:
        return returned_seqs
    else:
        return returned_seqs[0]
