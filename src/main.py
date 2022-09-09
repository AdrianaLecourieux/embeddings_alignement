# --- Main File --- #

import fonction  # import fonctions from the fonction file
import fasta_sequences 
import needlman_wunsch as NW
import sys # to interact with arguments

if __name__ == '__main__':
    
    seq1 = sys.argv[1] 
    seq2 = sys.argv[2]
    fasta1 = sys.argv[3]
    fasta2 = sys.argv[4]
    
    if not seq1.endswith(".t5emb") or not seq2.endswith(".t5emb"): 
        # Check the file extension
        raise Exception("ERROR : The two first arguments must be .t5emb extension")
    
    else: 
        embedding1_list = fonction.read_embedding(seq1)
        embedding2_list = fonction.read_embedding(seq2) 
        fasta1_list = fasta_sequences.read_fasta(fasta1)
        fasta2_list = fasta_sequences.read_fasta(fasta2)
                
        # Check the correspondance between embedding and fasta file 
        if len(embedding1_list) != len(fasta1_list) or len(embedding2_list) != len(fasta2_list):
            raise Exception("ERROR : For one sequence, embedding and fasta sequence must have the same length")

        # Calcul of the dot product matrix
        dot_pro_mat = fonction.dot_product(embedding1_list, embedding2_list)
        
        # Realise Needlman and Wunsh alignment 
        transformed_matrix, dot_matrix = NW.do(dot_pro_mat, fasta1_list, fasta2_list)
        
        output1, output2 = NW.needlman_wunsch(fasta1_list,fasta2_list, transformed_matrix, dot_matrix)
        print(output1 + "\n" + output2)
    


