##! /usr/bin/python3

# --- Pre-processing of fasta files

def read_fasta(fasta_file):
    """Open and read fasta file.

    Remove the description line.
    Separate each position with quotes.
    Add it to a list.

    Parameters
    ----------
    file : fasta file
       The first line is the protein description.
       The second line is the sequence of amino acids of the protein
               
    Returns
    -------
    list_fasta : list
        list with every amino acids of fasta_file separated by quotes 
    """
    line_sequence = ""
    with open(fasta_file, "r") as fasta:
        
        for line_in_fasta in fasta:
            
            if not line_in_fasta.startswith(">"):
                line_sequence += line_in_fasta.strip()
                
        list_fasta = list(line_in_fasta)
        
    return(list_fasta) 






