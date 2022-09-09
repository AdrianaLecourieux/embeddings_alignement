##! /usr/bin/python3
"""Function to put a sequence from a fasta file to list. (sorry pr les fautes)"""

def read_fasta(file):
    line_seq = ""
    with open(file, "r") as fasta:
        for line in fasta:
            
            if not line.startswith(">"):
                line_seq += line.strip()
                print(line_seq)
        seq = list(line_seq)
    return seq 


if __name__ == "__main__":
    fasta = "5_3_EXONUCLEASE_1BGXT.fasta"
    seq = read_fasta(fasta)
    print(len(seq))
    print(f'{seq}')

# du coup les fichiers ont généralement ce genre de format voilaa


 """   
    def check_len_fasta_emb (emb, fasta):
        if len(seq) != len(emb):
            print("ERROR : FASTA and Embeddings aren't corrresponding
                  "(not the same length)")
        else:
            return()
 """