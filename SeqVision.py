from abc import ABC, abstractmethod
from typing import Union


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
