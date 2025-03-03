import os
from typing import Union, Tuple
from abc import ABC, abstractmethod

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):
    @property
    @abstractmethod
    def seq(self) -> str:
        """An abstract property for obtaining a sequence"""
        pass

    @property
    @abstractmethod
    def allowed_chars(self) -> set:
        """An abstract property for getting valid letters"""
        pass

    def __len__(self) -> int:
        """Returns the length of the sequence"""
        return len(self.seq)

    def __getitem__(self, index: Union[int, slice]) -> Union[str, "BiologicalSequence"]:
        """Returns a slice of the sequence"""
        if isinstance(index, slice):
            sliced = self.seq[index]
            return self._create_new(sliced)
        return self.seq[index]

    def __str__(self) -> str:
        """Returns a string representation of the sequence"""
        return self.seq

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}: {self.seq}"

    def is_valid_alphabet(self) -> bool:
        """Сheck the alphabet of the sequence for correctness"""
        invalid_chars = {char for char in self.seq if char not in self.allowed_chars}
        if invalid_chars:
            raise ValueError(f"Invalid characters in the sequence: {invalid_chars}")
        return True

    @abstractmethod
    def _create_new(self, new_seq: str) -> "BiologicalSequence":
        pass


class NucleicAcidSequence(BiologicalSequence, ABC):
    @property
    @abstractmethod
    def _complement_dict(self) -> dict:
        pass

    def complement(self) -> "NucleicAcidSequence":
        return self._create_new(
            "".join(self._complement_dict[base] for base in self.seq)
        )

    def reverse(self) -> "NucleicAcidSequence":
        return self._create_new(self.seq[::-1])

    def reverse_complement(self) -> "NucleicAcidSequence":
        return self.reverse().complement()


class DNASequence(NucleicAcidSequence):
    def __init__(self, seq: str):
        self._seq = seq.upper()
        self.is_valid_alphabet()

    _complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    allowed_chars = {"A", "T", "C", "G"}

    @property
    def _transcribe_dict(self) -> dict:
        return {"A": "U", "T": "A", "C": "G", "G": "C"}

    @property
    def seq(self) -> str:
        return self._seq

    def transcribe(self) -> "RNASequence":
        return RNASequence(self._seq.translate(str.maketrans(self._transcribe_dict)))

    def _create_new(self, new_seq: str) -> "DNASequence":
        return DNASequence(new_seq)


class RNASequence(NucleicAcidSequence):
    def __init__(self, seq: str):
        self._seq = seq.upper()
        self.is_valid_alphabet()

    _complement_dict = {"A": "U", "U": "A", "C": "G", "G": "C"}
    allowed_chars = {"A", "U", "C", "G"}

    @property
    def seq(self) -> str:
        return self._seq

    def _create_new(self, new_seq: str) -> "RNASequence":
        return RNASequence(new_seq)


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, seq: str):
        self._seq = seq
        self.is_valid_alphabet()

    _MOLECULAR_WEIGHTS = {
        "A": 89.09,
        "R": 174.20,
        "N": 132.12,
        "D": 133.10,
        "C": 121.16,
        "Q": 146.15,
        "E": 147.13,
        "G": 75.07,
        "H": 155.16,
        "I": 131.18,
        "L": 131.18,
        "K": 146.19,
        "M": 149.21,
        "F": 165.19,
        "P": 115.13,
        "S": 105.09,
        "T": 119.12,
        "W": 204.23,
        "Y": 181.19,
        "V": 117.15,
    }
    allowed_chars = {
        "A",
        "R",
        "N",
        "D",
        "C",
        "Q",
        "E",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
    }

    @property
    def seq(self) -> str:
        return self._seq

    def get_molecular_weight(self) -> int:
        return sum(self._MOLECULAR_WEIGHTS[aa] for aa in self._seq)

    def _create_new(self, new_seq: str) -> "AminoAcidSequence":
        return AminoAcidSequence(new_seq)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[int, int], int, float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: Union[int, float] = 0,
) -> None:
    gc_min, gc_max = (
        (0, gc_bounds) if isinstance(gc_bounds, (int, float)) else gc_bounds
    )
    len_min, len_max = (
        (0, length_bounds) if isinstance(length_bounds, int) else length_bounds
    )

    input_path = os.path.join("data", input_fastq)
    output_path = os.path.join("filtered", output_fastq)

    with open(input_path) as in_handle, open(output_path, "w") as out_handle:
        filtered_records = (
            record
            for record in SeqIO.parse(in_handle, "fastq")
            if (
                len_min <= len(record.seq) <= len_max
                and (
                    lambda s: (
                        (gc_min <= (gc_fraction(s, ambiguous="ignore") * 100) <= gc_max)
                        if s
                        else False
                    )
                )(str(record.seq))
                and (lambda q: sum(q) / len(q) >= quality_threshold if q else False)(
                    record.letter_annotations.get("phred_quality", [])
                )
            )
        )

        SeqIO.write(filtered_records, out_handle, "fastq")
