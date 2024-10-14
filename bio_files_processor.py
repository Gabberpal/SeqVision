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
    output_fasta_path = os.path.join("filtered", output_fasta)
    with open(fasta_path, "r") as fasta_file:
        with open(output_fasta_path, "a") as output_fasta_file:
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

def parse_blast_output(input_blast: str, output_blast: str):
    """
    Parses a BLAST output file to extract and sort protein descriptions.

    Parameters:
    -----------
    input_blast : str
        Path to the input BLAST output file relative to the "data" directory.

    output_blast : str
        Path to the output file where the parsed and sorted protein descriptions will be written.

    Returns:
    --------
    None
        The function writes the parsed and sorted protein descriptions to the specified output file.

    Notes:
    ------
    - The function reads the input BLAST output file line by line.
    - It extracts protein descriptions starting from the line containing "Description".
    - The extracted descriptions are sorted alphabetically.
    - The sorted descriptions are written to the output file with each description followed by '...' and a newline.
    """
    blast_path = os.path.join("data", input_blast)
    output_blast_path = os.path.join("filtered", output_blast)
    with open(blast_path, "r") as blast_file:
        with open(output_blast_path, "a") as blast_output:
            flag = False
            result = []
            while True:
                line = blast_file.readline()
                if "Description" in line:
                    flag = True 
                    continue
                if flag:
                    if line.startswith('\n'):
                        break
                    result.append(line.split('...')[0])
            result.sort()
            for prot in result:
                blast_output.write(prot + '...' + '\n')

parse_blast_output("example_blast_results.txt", "blast_output.txt")
                 

