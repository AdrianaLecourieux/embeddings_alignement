# --- Saving alignments in .txt file

from pathlib import Path # to remove path

def save_in_txt(seq1, seq2, alignment_method, seq_aligned_list):
    """Save alignments in .txt file.

    Create and write a file, rename with : name of sequence 1, name of sequence2,
    method of alignment. When many alignments are optimal, there are separated with 
    two '\n' in the file

    Parameters
    ----------
    seq1: string
        name of the file given in argument
    seq2: string
        name of the file given in argument
    alignment_method : string
        name of the alignment method given in argument
    seq_aligned_list: list
        list of alignments
            
    Returns
    -------
    Returns nothing but print a message
    """    
    seq1_without_extension = Path(seq1).stem
    seq2_without_extension = Path(seq2).stem
    
    with open(f'../results/{seq1_without_extension}__{seq2_without_extension}_{alignment_method}_alignement.txt', "w") as file:
               
                for seq_aligned_1, seq_aligned_2 in seq_aligned_list:
                    
                    file.write(seq_aligned_1 + "\n" + seq_aligned_2 + "\n\n" )
    
    print("Alignment completed successfully !")
    
    
def save_in_txt_gp(seq1, seq2, alignment_method, seq_aligned_list):
    """Save alignments in .txt file.

    Create and write a file, rename with : name of sequence 1, name of sequence2,
    method of alignment. When many alignments are optimal, there are separated with 
    two '\n' in the file

    Parameters
    ----------
    seq1: string
        name of the file given in argument
    seq2: string
        name of the file given in argument
    alignment_method : string
        name of the alignment method given in argument
    seq_aligned_list: list
        list of alignments
               
    Returns
    -------
    Returns nothing but print a message
    """    
    seq1_without_extension = Path(seq1).stem
    seq2_without_extension = Path(seq2).stem
    
    with open(f'../results/{seq1_without_extension}__{seq2_without_extension}_{alignment_method}_affine_gap_alignement.txt', "w") as file:
               
                for seq_aligned_1, seq_aligned_2 in seq_aligned_list:
                    
                    file.write(seq_aligned_1 + "\n" + seq_aligned_2 + "\n\n" )
    
    print("Alignment completed successfully !")