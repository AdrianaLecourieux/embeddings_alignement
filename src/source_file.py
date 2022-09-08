# Source file


# Read the script arguments

import sys # to interact with arguments

def read_arg(seq1, seq2): 
    if seq1.endswith(".t5emb") or seq2.endswith(".t5emb") != True: 
        # Check the file extension
        print("ERROR : The two first arguments need to be .t5emb sequences")
    # if the extension is good, check if the dot product file exists
    else: 
    # call the check_dot_product function    
        check_dot_product(seq1, seq2)
        
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

if __name__ == '__main__':
    read_arg(sys.argv[1], sys.argv[2])


# Dot Product Calcul
"""

def is_number(n): #verifie si c'est un nombre
    try:
        float(n)   #essaye convertir
    except ValueError:
        return False
    return True


seq1 = open("5_3_exonuclease_1bgxt.t5emb", "r")
seq2 = open("6PF2K_1bif.t5emb", "r")
seq1splitline = seq1.read().splitlines() # recup ligne par ligne
seq2splitline = seq2.read().splitlines()
seq1.close()
seq2.close()
dotproductmatrix = []
for lineseq1 in seq1splitline: 
    dotproductline =[]
    for lineseq2 in seq2splitline:
        listseq = []
        listseq2 = []
        tmp = []
        for x in lineseq1.split(' '): # découpe ligne sur les espaces #NE PAS PRECISER
            if is_number(x): # (le découpage est pas considéré comme float) check que c'est un nombre
                tmp.append(float(x))
                listseq.append(tmp)
        tmp2 = []
        for x in lineseq2.split(' '):
            if is_number(x): 
                tmp2.append(float(x))
                listseq2.append(tmp2)
        listetranspose = np.transpose(listseq2)
        dotproduct = np.dot(listseq, listetranspose)
        dotproductline.append(dotproduct[0][0]) # récup que la valeur
    dotproductmatrix.append(dotproductline)
print(dotproductmatrix)

"""

