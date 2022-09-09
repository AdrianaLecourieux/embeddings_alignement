# --- Needlman and Wunsch alignment script

import numpy as np



seq1 = "AVR"
seq2 = "AR"



        


def do(dot_matrix, seq1, seq2):
    # Creation of the transformed matrix
    col_seq1 = len(seq1) + 1
    row_seq2 = len(seq2) + 1
    transformed_matrix= np.zeros((col_seq1, row_seq2),dtype = int)
    for i in range(1,col_seq1):
        
            for j in range(1,row_seq2):
                
                top = transformed_matrix[i][j-1] 
                diagonal = transformed_matrix[i-1][j-1] + dot_matrix[i-1][j-1] 
                bottom = transformed_matrix[i-1][j] 
                transformed_matrix[i,j] = max(diagonal, top, bottom)
                
    return(np.matrix(transformed_matrix), dot_matrix) 

def needlman_wunsch(seq1,seq2, transformed_matrix, dot_matrix):
    # Find optimal way
    aligned_sequence1 = ""
    aligned_sequence2 = ""
    print(transformed_matrix.shape)
    i = len(seq1)
    j = len(seq2)
    print(i, j)

    #From the right bottom to the left top
    while i > 0 and j > 0:
        print(i, j)     
        score_position = transformed_matrix[i,j]
        score_diagonal = transformed_matrix[i-1,j-1]
        score_top = transformed_matrix[i,j-1]
        score_left = transformed_matrix[i-1,j]
        
        tmp, tmpi, tmpj = score_diagonal, i-1 , j-1
        
        # Calcule the score value
        
        if score_position == score_top :
            aligned_sequence1 += seq1[j-1]
            aligned_sequence2 += '-'
            j -= 1
            
        elif score_position == score_diagonal + dot_matrix[i-1][j-1] :
            aligned_sequence1 += seq1[j-1]
            aligned_sequence2 += seq2[i-1]
            i -= 1
            j -= 1
                       
        elif score_position == score_left :
            print(i)
            aligned_sequence1 += '-'
            aligned_sequence2 += seq2[i-1]
            i -= 1
    # finish
    
    while j > 0:
        aligned_sequence1 += seq1[j-1]
        aligned_sequence2 += '-'
        j -= 1
        
    while i > 0:
        aligned_sequence1 += '-'
        aligned_sequence2 += seq2[i-1]
        i -= 1
    
    aligned_sequence1 = aligned_sequence1[::-1]
    aligned_sequence2 = aligned_sequence2[::-1]
    
    return(aligned_sequence1, aligned_sequence2)


    
