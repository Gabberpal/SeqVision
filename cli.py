import argparse
from bio_files_processor import convert_multiline_fasta_to_oneline, parse_blast_output
from SeqVision import DNASequence, RNASequence, AminoAcidSequence, filter_fastq


def main():
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
    filter_parser.add_argument("--length", type=int, nargs=2, metavar=("MIN", "MAX"))
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
    rna_parser.add_argument("--operation", choices=["complement", "reverse", "revcomp"])

    # Protein tools
    protein_parser = subparsers.add_parser("protein", help="Protein analysis")
    protein_parser.add_argument("--seq", required=True)

    args = parser.parse_args()

    if args.command == "filter_fastq":
        filter_fastq(
            input_fastq=args.input,
            output_fastq=args.output,
            gc_bounds=tuple(args.gc) if args.gc else (0, 100),
            length_bounds=tuple(args.length) if args.lenth else (0, 2**32),
            quality_threshold=args.quality,
        )

    elif args.command == "convert_fasta":
        convert_multiline_fasta_to_oneline(args.input, args.output)

    elif args.command == "parse_blast":
        parse_blast_output(args.input, args.output)

    elif args.command == "dna":
        dna = DNASequence(args.seq)
        if args.operation == "complement":
            print(dna.complement())
        elif args.operation == "reverse":
            print(dna.reverse())
        elif args.operation == "revcomp":
            print(dna.reverse_complement())
        elif args.operation == "transcribe":
            print(dna.transcribe())
        elif args.operation == "slice":
            print(dna[args.slice[0] : args.slice[1]])

    elif args.command == "rna":
        rna = RNASequence(args.seq)
        if args.operation == "complement":
            print(rna.complement())
        elif args.operation == "reverse":
            print(rna.reverse())
        elif args.operation == "revcomp":
            print(rna.reverse_complement())

    elif args.command == "protein":
        protein = AminoAcidSequence(args.seq)
        print(protein.get_molecular_weight())


if __name__ == "__main__":
    main()
