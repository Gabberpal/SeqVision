set_rna = {"A", "G", "C", "U", "a", "g", "c", "u"}
set_dna = {"A", "G", "C", "T", "a", "g", "c", "t"}
ComplementDictDna = {
    "a": "t",
    "A": "T",
    "t": "a",
    "T": "A",
    "g": "c",
    "G": "C",
    "c": "g",
    "C": "G",
}
ComplementDictRna = {
    "a": "u",
    "A": "U",
    "u": "a",
    "U": "A",
    "g": "c",
    "G": "C",
    "c": "g",
    "C": "G",
}
TranscribeDictDna = {
    "a": "a",
    "A": "A",
    "t": "u",
    "T": "U",
    "g": "g",
    "G": "G",
    "c": "c",
    "C": "C",
}
TranscribeDictRna = {
    "a": "a",
    "A": "A",
    "u": "t",
    "U": "T",
    "g": "g",
    "G": "G",
    "c": "c",
    "C": "C",
}


def is_rna(seq):
    if set(seq.upper()) <= set_rna:
        return True
    else:
        return False


def is_dna(seq):
    if set(seq.upper()) <= set_dna:
        return True
    else:
        return False


def complement(seq):
    complement_seq = ""
    if is_dna(seq):
        complement_seq = "".join([ComplementDictDna[i] for i in seq])
    if is_rna(seq):
        for i in seq:
            complement_seq = "".join([ComplementDictRna[i] for i in seq])
    return complement_seq


def reverse(seq):
    return seq[::-1]


def transcribe(seq):
    transcribe_seq = ""
    if is_dna(seq):
        transcribe_seq = "".join([TranscribeDictDna[i] for i in seq])
    if is_rna(seq):
        transcribe_seq = "".join([TranscribeDictRna[i] for i in seq])
    return transcribe_seq


def reverse_complement(seq):
    return reverse(complement(seq))


# Функция, которая возвращает последовательность, если она является палиндромом.
# Удобно, если из большого числа последовательностей нужно отобрать только палиндромы
def which_palindrome(seq):
    if len(seq) % 2 == 0:
        if (
            seq[: int(len(seq) / 2)].lower()
            == seq[-1 : int(len(seq) / 2 - 1) : -1].lower()
        ):
            return seq


functions = {
    "transcribe": transcribe,
    "complement": complement,
    "reverse": reverse,
    "reverse_complement": reverse_complement,
    "which_palindrome": which_palindrome,
}
