# --- Main File --- #

from ctypes import alignment
import fonction  # import fonctions from the fonction file
import fasta_sequences 
import needlman_wunsch as NW
import smith_waterman as SW
import os
import sys
import argparse  # to interact with arguments


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    
    # Required arguments
    
    parser.add_argument("-emb1","--embedding1",  help= "Embedding 1 in .t5emb extension")
    parser.add_argument("-emb2","--embedding2",  help= "Embedding 2 in .t5emb extension")
    parser.add_argument("-f1","--fasta1",  help= "Fasta 1 in .FASTA extension")
    parser.add_argument("-f2","--fasta2",  help= "Fasta 2 in .FASTA extension")
    
    # Optionnal argument
    
    parser.add_argument("-m","--method",  help='Choose a "global" (Needleman and Wunsch) or "local" (Smith and Waterman) alignment algorithm.  -m global default"')
    
    # Assign arguments to variables
    
    args = parser.parse_args()    
    seq1 = args.embedding1
    seq2 = args.embedding2
    fasta1 = args.fasta1
    fasta2 = args.fasta2
    alignment_method = args.method

    # Check the file extension
     
    if not seq1.endswith(".t5emb") or not seq2.endswith(".t5emb"):   
        raise Exception("ERROR : The two first arguments must be .t5emb extension")
    
    
    else: 
        embedding1_list = fonction.read_embedding(seq1)
        embedding2_list = fonction.read_embedding(seq2) 
        fasta1_list = fasta_sequences.read_fasta(fasta1)
        fasta2_list = fasta_sequences.read_fasta(fasta2)
        
                
        # Check the correspondance between embedding and fasta file 
        
        if len(embedding1_list) != len(fasta1_list) or len(embedding2_list) != len(fasta2_list):
            raise Exception("ERROR : For one sequence, embedding and fasta sequence must have the same length.")

        # Calcul of the dot product matrix
        
        dot_pro_mat = fonction.dot_product(embedding1_list, embedding2_list)
        
        
        if alignment_method == "global" or not alignment_method :
            
            # ------ Realise Needlman and Wunsh alignment 
            
            transformed_matrix = NW.transformation_NW(dot_pro_mat, fasta1_list, fasta2_list)
            
            seq_aligned_1, seq_aligned_2 = NW.needlman_wunsch(fasta1_list,fasta2_list, transformed_matrix)
            print("\n" + seq_aligned_1 + "\n" + seq_aligned_2 + "\n" + "Alignment completed successfully !" )

            # Save output in .txt file
            
            seq1_without_extenstion = os.path.splitext(seq1)[0]
            seq2_without_extension = os.path.splitext(seq2)[0]
            
            with open(f'../results/{seq1_without_extenstion}__{seq2_without_extension}_global_alignement.txt', "w") as file:
                file.write(seq_aligned_1 + "\n" + seq_aligned_2)
        
        
        elif alignment_method == "local":
            
            # ------ Realise Smith and Waterman alignment 
            
            transformed_matrix_SW = SW.transformation_SW(dot_pro_mat, fasta1_list, fasta2_list)
            seq_aligned_1_SW, seq_aligned_2_SW = SW.smith_waterman(fasta1_list,fasta2_list, transformed_matrix_SW)
            print("\n" + seq_aligned_1_SW + "\n" + seq_aligned_2_SW + "\n" + "Alignment completed successfully !" )
            
            # Save output in .txt file
            
            seq1_without_extenstion_SW = os.path.splitext(seq1)[0]
            seq2_without_extension_SW = os.path.splitext(seq2)[0]
            with open(f'../results/{seq1_without_extenstion_SW}__{seq2_without_extension_SW}_local_alignement.txt', "w") as file:
                file.write(seq_aligned_1_SW + "\n" + seq_aligned_2_SW)
        