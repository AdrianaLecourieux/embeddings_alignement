import numpy as np

def transformation_semi_global(dot_matrix, seq1, seq2):
    
    # Creation of the transformed matrix
    row_seq1 = len(seq1) + 1
    col_seq2 = len(seq2) + 1
    
    transformed_matrix= np.zeros((row_seq1, col_seq2),dtype = int)
    
    # Case where we want penalty-free in the first col, so fix the first col to 0
    transformed_matrix[:, 0] =  0
    
    # Case where we want penalty-free in the first row, so fix the first row to 0
    #transformed_matrix[0,:] =  0
    
    
    for i in range(0,row_seq1):
        
            for j in range(1,col_seq2):
                
                top = transformed_matrix[i][j-1] 
                diagonal = transformed_matrix[i-1][j-1] + dot_matrix[i-1][j-1] 
                left = transformed_matrix[i-1][j] 
                transformed_matrix[i,j] = max(diagonal, top, left)
                
    return(np.matrix(transformed_matrix)) 



def semi_global(seq1,seq2, transformed_matrix):

    # Find optimal path
    
    aligned_sequence1 = ""
    aligned_sequence2 = ""
    
    # we want the begin with the best score in the last column # AAAAAAAAAAAHH
    # the max of i in j col
    i = len(seq1)
    j = len(seq2)
    
    j = np.amax(transformed_matrix[:,j], axis=0)

    #From the last column to the first colum,
    
    while i > 0 and j > 0:
        
        score_diagonal = transformed_matrix[i-1,j-1]
        score_top = transformed_matrix[i,j-1]
        score_left = transformed_matrix[i-1,j]
        max_value = max(score_diagonal, score_left, score_top)        

        # Calcule the max score value
            
        if max_value == score_diagonal :

            aligned_sequence2 += seq2[j-1]
            aligned_sequence1 += seq1[i-1]
            i -= 1
            j -= 1
        
        elif max_value == score_top :

            aligned_sequence2 += seq2[j-1]
            aligned_sequence1 += '-'
            j -= 1
                       
        elif max_value == score_left :
            
            aligned_sequence2 += '-'
            aligned_sequence1 += seq1[i-1]
            i -= 1
            
    # finish when meeting 0 on j or i
    
    while j > 0:
        
        aligned_sequence2 += seq2[j-1]
        aligned_sequence1 += '-'
        j -= 1
        
    while i > 0:
        
        aligned_sequence2 += '-'
        aligned_sequence1 += seq1[i-1]
        i -= 1
    
    aligned_sequence1 = aligned_sequence1[::-1]
    aligned_sequence2 = aligned_sequence2[::-1]
    
    return(aligned_sequence1, aligned_sequence2)

