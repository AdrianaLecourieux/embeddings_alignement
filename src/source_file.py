# Source file

"""
with open("../data/fasta_sequences/5_3_EXONUCLEASE_1BGXT.fasta") as f:
    firstline = f.readline().rstrip()

print(firstline)
"""

"""
l1 = [1,2,3]
l2 = [4,5,6]
l3 = [7,8,9]

seq1 = []
seq1.append(l1)
seq1.append(l2)
seq1.append(l3)



v1 = [1,2,3]
v2 = [4,5,6]
v3 = [7,8,9]

seq2 = []
seq2.append(v1)
seq2.append(v2)
seq2.append(v3)

print(seq1, seq2)
"""
import numpy as np

seq1 = [1,2,3]
seq2 = [1,2,3]

dotproduct = np.dot(seq1, seq2)

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
        for x in lineseq1.split(' '): # découpe ligne sur les espaces
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

listseq = []
for y in lignes1:
    tmp = []
    for x in y.split(' '):
        tmp.append(int(x))
    listseq.append(tmp)
print(listseq)




lignes2 = seq2.read().splitlines()
listseq2 = []
for y in lignes2:
    tmp2 = []
    for x in y.split(' '):
        tmp2.append(int(x))
    listseq2.append(tmp2)
print(listseq2)

listetranspose = np.transpose(listseq2)

dotproduct = np.dot(listseq, listetranspose)
print(dotproduct)

"""