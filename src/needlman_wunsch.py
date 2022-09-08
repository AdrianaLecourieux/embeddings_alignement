# --- Needlman and Wunsch alignment script

import numpy as np

from fonction import dot_product
"""
seq1 = "ATGC"
seq2 = "GGATGC"


length_seq1 = len(seq1)
length_seq2 = len(seq2)
gap = 0
F = np.zeros((length_seq1 + 1, length_seq2+ 1))
print(F)

F[:,0] = np.linspace(0, -nx * gap, nx + 1)
print(F)

def nw(x, y, match = 1, mismatch = 1, gap = 1):
    nx = len(x)
    ny = len(y)
    # Optimal score at each possible pair of characters.
    F = np.zeros((nx + 1, ny + 1))
"""
seq1 = "AVR"
seq2 = "AR"

col_seq1 = len(seq1) + 1
row_seq2 = len(seq2) + 1

        


dot_matrix = [[1, 2], [1, 2], [1, 3]]

def do(dot_matrix):
    transformed_matrix= np.zeros((col_seq1, row_seq2),dtype = int)
    for i in range(1,col_seq1):
            for j in range(1,row_seq2):
                diagonal = transformed_matrix[i-1][j-1] + dot_matrix[i-1][j-1] 
                top = transformed_matrix[i][j-1] 
                bottom = transformed_matrix[i-1][j] 
                transformed_matrix[i,j] = max(diagonal, top, bottom)
    return(np.matrix(transformed_matrix)) 



