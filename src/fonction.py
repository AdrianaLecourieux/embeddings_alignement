# --- Fonction file

# This script contains the functions used from the beginning of the program to 
# the alignement

## --- Read the script arguments

import sys
from turtle import shape # to interact with arguments
import numpy as np

def read_arg(seq1, seq2): 
    
    if not seq1.endswith(".t5emb") or not seq2.endswith(".t5emb"): 
        # Check the file extension
        print("ERROR : The two first arguments need to be .t5emb extension")
    
    else: 
        embedding1_list = read_embedding(seq1)
        embedding2_list = read_embedding(seq2) 
        dot_product(embedding1_list, embedding2_list)
    return(embedding1_list, embedding2_list)
    
 
      
def read_embedding(file):
    embedding_list = []
    with open(file, "r") as embedding:
        
        for line in embedding:
            
            vector = line.split()
            vector = [float(x) for x in vector]
            embedding_list.append(vector)  
               
    return(embedding_list)

#comment utiliser la fonction

def dot_product(emb1_list, emb2_list):
    calcul_dot_product = np.dot(emb1_list, np.array(emb2_list).T)
    return(calcul_dot_product)
   # np.savetxt(file_out, res, delimiter="\t")

    
