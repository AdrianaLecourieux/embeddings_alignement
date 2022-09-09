# --- Main File --- #

import fonction as fonction # import fonctions from the fonction file
import needlman_wunsch as NW
import sys # to interact with arguments

if __name__ == '__main__':
    
    seq1 = sys.argv[1] 
    seq2 = sys.argv[2]
    fasta1 = sys.argv[3]
    fasta2 = sys.argv[4]
    
    if not seq1.endswith(".t5emb") or not seq2.endswith(".t5emb"): 
        # Check the file extension
        print("ERROR : The two first arguments need to be .t5emb extension")

    
    else: 
        embedding1_list = fonction.read_embedding(seq1)
        embedding2_list = fonction.read_embedding(seq2) 
        dot_pro_mat = fonction.dot_product(embedding1_list, embedding2_list)
        
        alignment = NW.do(dot_pro_mat)
        print(alignment)
    


