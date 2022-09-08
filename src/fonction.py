# --- Fonction file

# This script contains the functions used from the beginning of the program to 
# the alignement

## Read the script arguments

import sys # to interact with arguments

def read_arg(seq1, seq2): 
    
    if seq1.endswith(".t5emb") or seq2.endswith(".t5emb") != True: 
        # Check the file extension
        print("ERROR : The two first arguments need to be .t5emb sequences")
    # if the extension is good, check if the dot product file exists
    else: 
        
    # call the check_dot_product function    
        check_dot_product(seq1, seq2)
      



## Check if the dot product matrix file exists

import os.path # to interact with the operating system

def check_dot_product(emb1, emb2):
    emb1split = os.path.splitext(emb1)
    emb2split = os.path.splitext(emb2)
    path_to_doproduct = f'../results/{emb1split[0]}__{emb2split[0]}.txt'
    
    if os.path.exists(path_to_doproduct):
        print("The dot product matrice exists")
        go_to_alignement_matrice(emb1, emb2)
        
    else:
        print("The dot product matrice doesn't exist") 
        go_to_dot_product(emb1, emb2) # 
        
           
def go_to_dot_product(input1, input2):
    print(input1, input2)

def go_to_alignement_matrice(input1, input2):
    print("coucou")