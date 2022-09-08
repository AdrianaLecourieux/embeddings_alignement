##! /usr/bin/python3
"""Function to put a sequence from a fasta file to list. (sorry pr les fautes)"""

def read_fasta(file):
    """Get a fasta file as entry to return the corresponding seq in a list.

    Parameters
    ----------
    file : str 
        Name of fasta file.

    Returns
    -------
    list
        Sequence with each position separated.
    """
    line_seq = ""
    with open(file, "r") as fasta:
        for line in fasta:
            if not line.startswith(">"):
                line_seq += line.strip()
                print(line_seq)
        seq = list(line_seq)
    return seq 


if __name__ == "__main__":
    fasta = "P00509.fasta"
    seq = read_fasta(fasta)
    print(f'{seq}')

# du coup les fichiers ont généralement ce genre de format voilaa