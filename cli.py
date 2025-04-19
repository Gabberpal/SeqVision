import argparse
import logging
from bio_files_processor import convert_multiline_fasta_to_oneline, parse_blast_output
from SeqVision import DNASequence, RNASequence, AminoAcidSequence, filter_fastq

# Logger config
logging.basicConfig(
    filename="seqvision.log",
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger()


def main():
    try:
        logger.info("Запуск программы SeqVision...")

        parser = argparse.ArgumentParser(
            prog="SeqVision", description="Biological sequence toolkit"
        )
        subparsers = parser.add_subparsers(dest="command", required=True)

        # Filter FASTQ
        filter_parser = subparsers.add_parser(
            "filter_fastq", help="Filter FASTQ by GC, length and quality"
        )
        filter_parser.add_argument("--input", required=True)
        filter_parser.add_argument("--output", required=True)
        filter_parser.add_argument("--gc", type=float, nargs=2, metavar=("MIN", "MAX"))
        filter_parser.add_argument(
            "--length", type=int, nargs=2, metavar=("MIN", "MAX")
        )
        filter_parser.add_argument("--quality", type=float, default=0)

        # Conver FASTA
        fasta_parser = subparsers.add_parser(
            "convert_fasta", help="Convert multiline FASTA to single-line"
        )
        fasta_parser.add_argument("--input", required=True)
        fasta_parser.add_argument("--output", required=True)

        # Parse BLAST
        blast_parser = subparsers.add_parser(
            "parse_blast", help="Parse BLAST output and sort descriptions"
        )
        blast_parser.add_argument("--input", required=True)
        blast_parser.add_argument("--output", required=True)

        # DNA Tools
        dna_parser = subparsers.add_parser("dna", help="DNA sequence operations")
        dna_parser.add_argument("--seq", required=True)
        dna_parser.add_argument(
            "--operation",
            choices=["complement", "reverse", "revcomp", "transcribe", "slice"],
        )
        dna_parser.add_argument("--slice", type=int, nargs=2, metavar=("START", "END"))

        # RNA Tools
        rna_parser = subparsers.add_parser("rna", help="RNA sequence operations")
        rna_parser.add_argument("--seq", required=True)
        rna_parser.add_argument(
            "--operation", choices=["complement", "reverse", "revcomp"]
        )

        # Protein tools
        protein_parser = subparsers.add_parser("protein", help="Protein analysis")
        protein_parser.add_argument("--seq", required=True)

        args = parser.parse_args()

        if args.command == "filter_fastq":
            logger.info(f"Фильтрация FASTQ: {args.input} -> {args.output}")
            filter_fastq(
                input_fastq=args.input,
                output_fastq=args.output,
                gc_bounds=tuple(args.gc) if args.gc else (0, 100),
                length_bounds=tuple(args.length) if args.lenth else (0, 2**32),
                quality_threshold=args.quality,
            )
            logger.info(f"Фильтрация завершена для файла {args.input}")

        elif args.command == "convert_fasta":
            logger.info(f"Конвертация FASTA: {args.input} -> {args.output}")
            convert_multiline_fasta_to_oneline(args.input, args.output)
            logger.info(f"Конвертация завершена для файла {args.input}")

        elif args.command == "parse_blast":
            logger.info(f"Обработка BLAST результатов: {args.input} -> {args.output}")
            parse_blast_output(args.input, args.output)
            logger.info(f"Обработка BLAST завершена для файла {args.input}")

        elif args.command == "dna":
            dna = DNASequence(args.seq)
            logger.info(f"Операция с ДНК: {args.seq} ({args.operation})")
            if args.operation == "complement":
                result = dna.complement()
            elif args.operation == "reverse":
                result = dna.reverse()
            elif args.operation == "revcomp":
                result = dna.reverse_complement()
            elif args.operation == "transcribe":
                result = dna.transcribe()
            elif args.operation == "slice":
                result = dna[args.slice[0] : args.slice[1]]
            logger.info(f"Результат операции с ДНК: {result}")
            print(result)

        elif args.command == "rna":
            rna = RNASequence(args.seq)
            logger.info(f"Операция с РНК: {args.seq} ({args.operation})")
            if args.operation == "complement":
                result = rna.complement()
            elif args.operation == "reverse":
                result = rna.reverse()
            elif args.operation == "revcomp":
                result = rna.reverse_complement()
            elif args.operation == "transcribe":
                result = rna.transcribe()
            elif args.operation == "slice":
                result = rna[args.slice[0] : args.slice[1]]
            logger.info(f"Результат операции с РНК: {result}")
            print(result)

        elif args.command == "protein":
            protein = AminoAcidSequence(args.seq)
            logger.info(f"Операция с белком: {args.seq}")
            result = protein.get_molecular_weight()
            logger.info(f"Молекулярная масса белка: {result}")
            print(result)

    except Exception as e:
        logger.error(f"Ошибка: {e}")
        raise


if __name__ == "__main__":
    main()
