import os


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str):
    """
    Converts a multiline FASTA file to a single-line FASTA file.

    Parameters:
    -----------
    input_fasta : str
        Path to the input multiline FASTA file relative to the "data" directory.

    output_fasta : str
        Path to the output single-line FASTA file.

    Returns:
    --------
    None
        The function writes the converted FASTA data to the specified output file.

    Notes:
    ------
    - The function reads the input FASTA file line by line.
    - It writes each sequence on a single line in the output file.
    - The function handles the case where the input file contains sequences that span multiple lines.
    """
    fasta_path = os.path.join("data", input_fasta)
    with open(fasta_path, "r") as fasta_file:
        with open(output_fasta, "a") as output_fasta_file:
            prev = "0"
            while True:
                line = fasta_file.readline().strip()
                if not line:
                    break
                elif line.startswith(">") and prev == "0":
                    output_fasta_file.write(line + "\n")
                    prev = line
                    continue
                elif (not prev.startswith(">")) and line.startswith(">"):
                    output_fasta_file.write("\n" + line + "\n")
                    prev = line
                    continue
                else:
                    output_fasta_file.write(line)
                    prev = line


