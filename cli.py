import argparse
import logging
from bio_files_processor import convert_multiline_fasta_to_oneline, parse_blast_output
from SeqVision import DNASequence, RNASequence, AminoAcidSequence, filter_fastq

# Logger configuration
logging.basicConfig(
    filename="seqvision.log",
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger()


def main():
    try:
        logger.info("Starting SeqVision CLI...")

        parser = argparse.ArgumentParser(
            prog="SeqVision", description="Biological sequence toolkit"
        )
        subparsers = parser.add_subparsers(dest="command", required=True)

        # FASTQ Filtering
        filter_parser = subparsers.add_parser(
            "filter_fastq", help="Filter FASTQ by GC content, length, and quality"
        )
        filter_parser.add_argument("--input", required=True)
        filter_parser.add_argument("--output", required=True)
        filter_parser.add_argument("--gc", type=float, nargs=2, metavar=("MIN", "MAX"))
        filter_parser.add_argument(
            "--length", type=int, nargs=2, metavar=("MIN", "MAX")
        )
        filter_parser.add_argument("--quality", type=float, default=0)

        # FASTA Conversion
        fasta_parser = subparsers.add_parser(
            "convert_fasta", help="Convert multiline FASTA to single-line format"
        )
        fasta_parser.add_argument("--input", required=True)
        fasta_parser.add_argument("--output", required=True)

        # BLAST Parsing
        blast_parser = subparsers.add_parser(
            "parse_blast", help="Parse BLAST output and sort protein descriptions"
        )
        blast_parser.add_argument("--input", required=True)
        blast_parser.add_argument("--output", required=True)

        # DNA Operations
        dna_parser = subparsers.add_parser(
            "dna", help="Perform operations on DNA sequences"
        )
        dna_parser.add_argument("--seq", required=True)
        dna_parser.add_argument(
            "--operation",
            choices=["complement", "reverse", "revcomp", "transcribe", "slice"],
        )
        dna_parser.add_argument("--slice", type=int, nargs=2, metavar=("START", "END"))

        # RNA Operations
        rna_parser = subparsers.add_parser(
            "rna", help="Perform operations on RNA sequences"
        )
        rna_parser.add_argument("--seq", required=True)
        rna_parser.add_argument(
            "--operation", choices=["complement", "reverse", "revcomp"]
        )

        # Protein Analysis
        protein_parser = subparsers.add_parser(
            "protein", help="Analyze protein sequence"
        )
        protein_parser.add_argument("--seq", required=True)

        args = parser.parse_args()

        if args.command == "filter_fastq":
            logger.info(f"Filtering FASTQ: {args.input} -> {args.output}")
            filter_fastq(
                input_fastq=args.input,
                output_fastq=args.output,
                gc_bounds=tuple(args.gc) if args.gc else (0, 100),
                length_bounds=tuple(args.length) if args.length else (0, 2**32),
                quality_threshold=args.quality,
            )
            logger.info(f"Filtering completed for file {args.input}")

        elif args.command == "convert_fasta":
            logger.info(f"Converting FASTA: {args.input} -> {args.output}")
            convert_multiline_fasta_to_oneline(args.input, args.output)
            logger.info(f"FASTA conversion completed for file {args.input}")

        elif args.command == "parse_blast":
            logger.info(f"Parsing BLAST output: {args.input} -> {args.output}")
            parse_blast_output(args.input, args.output)
            logger.info(f"BLAST output processing completed for file {args.input}")

        elif args.command == "dna":
            dna = DNASequence(args.seq)
            logger.info(f"DNA operation: {args.seq} ({args.operation})")
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
            logger.info(f"DNA operation result: {result}")
            print(result)

        elif args.command == "rna":
            rna = RNASequence(args.seq)
            logger.info(f"RNA operation: {args.seq} ({args.operation})")
            if args.operation == "complement":
                result = rna.complement()
            elif args.operation == "reverse":
                result = rna.reverse()
            elif args.operation == "revcomp":
                result = rna.reverse_complement()
            logger.info(f"RNA operation result: {result}")
            print(result)

        elif args.command == "protein":
            protein = AminoAcidSequence(args.seq)
            logger.info(f"Protein analysis: {args.seq}")
            result = protein.get_molecular_weight()
            logger.info(f"Protein molecular weight: {result}")
            print(result)

    except Exception as e:
        logger.error(f"An error occurred: {e}")
        raise


if __name__ == "__main__":
    main()
