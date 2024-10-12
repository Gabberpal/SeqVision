# SeqVision 

**SeqVision** - это инструмент для обработки и анализа нуклеотидных последовательностей, таких как ДНК и РНК. Он предоставляет функции для фильтрации, трансформации и анализа последовательностей.

## Оглавление 

- [Функциональность](#функциональность)
- [Установка](#установка)
- [Использование](#использование)
- [Примеры](#примеры)

## Функциональность 

- **Фильтрация последовательностей**: Фильтрация последовательностей по GC-составу, длине и качеству.
- **Трансформация последовательностей**: Транскрипция, комплементарность, обратная комплементарность, реверс и отбор палиндромных последовательностей.
- **Анализ последовательностей**: Вычисление статистики по последовательностям, такой как GC-состав и длина.

## Установка 

Клонируйте репозиторий:

   ```bash
   git clone https://github.com/Gabberpal/SeqVision
   cd SeqVision
   ```

## Использование 

### Фильтрация последовательностей

Для фильтрации последовательностей по GC-составу, длине и качеству, используйте функцию `filter_fastq`. Пример:

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
### Трансформация последовательностей 

```python
sequences = ["ACGT", "UGCA"]
transformed_sequences = run_dna_rna_tools(sequences, "function_name")
print(transformed_sequences)
```

#### Функции для трансформации последовательностей 

- `transcribe`: возвращает транскрибированные последовательности.
- `complement`: возвращает комплементарные последовательности.
- `reverse`: возвращает обратные последовательности.
- `reverse_complement`: возвращает обратно-комплементарные последовательности.
- `which_palindrome`: возврощает только палиндромные последовательности.


