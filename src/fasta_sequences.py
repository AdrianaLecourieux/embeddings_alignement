##! /usr/bin/python3


def read_fasta(fasta_file):

    line_sequence = ""
    with open(fasta_file, "r") as fasta:
        
        for line_in_fasta in fasta:
            
            if not line_in_fasta.startswith(">"):
                line_sequence += line_in_fasta.strip()
                
        list_fasta = list(line_in_fasta)
        #print(f'{list_fasta}')
        
    return(list_fasta) 



# du coup les fichiers ont généralement ce genre de format voilaa


