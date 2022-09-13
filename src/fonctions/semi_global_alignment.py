# --- Needleman and Wunsch alignment script

from re import M
import numpy as np # import the numpy module and rename it

def transformation_semi_global(dot_matrix, seq1, seq2):
    """Creation of the transformed matrix from the dot product matrix and 
    fasta sequences.

    The transformed matrix is created by following the semi-global algorithm (based
    on NW algorithm) 
    with the dot product matrix as the score matrix and a fixed gap penalty set to 0

    Parameters
    ----------
    dot_matrix : array
        An array where each value corresponds to the dot product obtained 
        from the embedding vectors of one position of protein 1 and another
        position of protein 2.
    seq2: list
        list with every amino acids of fasta_file2 separated by quotes
    seq2: list
        list with every amino acids of fasta_file2 separated by quotes
        
    Returns
    -------
    transformed_matrix : matrix
        New matrix completed from the semi-global algorithm 
    """
    # Creation of the transformed matrix
    row_seq= len(seq2) + 1
    col_seq = len(seq1) + 1
    
    transformed_matrix= np.zeros((row_seq, col_seq),dtype = int)
    
    # Case where we want penalty-free in the first col, so fix the first col to 0
    transformed_matrix[:, 0] =  0
    
    for i in range(0,row_seq):
        
            for j in range(1,col_seq):
                
                left = transformed_matrix[i][j-1] 
                diagonal = transformed_matrix[i-1][j-1] + dot_matrix[i-1][j-1] 
                top = transformed_matrix[i-1][j] 
                transformed_matrix[i,j] = max(diagonal, top, left)
            
    return(np.matrix(transformed_matrix)) 



def semi_global(seq1,seq2, transformed_matrix):
    """Semi-global alignment.

    Finds the optimal path starting from maximum value in the last column of the transformed 
    matrix and going up to the first column . Creating chain of string with gap when insertion
    or deletion occurs (top / left).

    Parameters
    ----------
    seq2: list
        list with every amino acids of fasta_file2 separated by quotes
    seq2: list
        list with every amino acids of fasta_file2 separated by quotes
    transformed_matrix : matrix
        Transformed matrix completed from the Needleman and Wunsch algorithm 
          
    Returns
    -------
    aligned_sequence1: string
        Chain of string composed of amino acid and gap
    aligned_sequence2: string
        Chain of string composed of amino acid and gap
    """
    # Find optimal path
    
    aligned_sequence1 = ""
    aligned_sequence2 = ""
    
    # we want the begin with the best score in the last column 
    # the max of i in j col
    j = len(seq1)
    index_max = np.where(transformed_matrix[:,j] == np.amax(transformed_matrix[:,j]))
    
    #i = int(index_max[0])
    idx_row = index_max[0]
    i = int(idx_row[0])
    print(i)
    #From the last column to the first colum,
    
    while i > 0 and j > 0:
        
        score_diagonal = transformed_matrix[i-1,j-1]
        score_left = transformed_matrix[i,j-1]
        score_top = transformed_matrix[i-1,j]
        max_value = max(score_diagonal, score_left, score_top)        

        # Calcule the max score value
            
        if max_value == score_diagonal :

            aligned_sequence2 += seq2[i-1]
            aligned_sequence1 += seq1[j-1]
            i -= 1
            j -= 1
        
        elif max_value == score_top:

            aligned_sequence2 += seq2[i-1]
            aligned_sequence1 += '-'
            i -= 1
                       
        elif max_value == score_left :
            
            aligned_sequence2 += '-'
            aligned_sequence1 += seq1[j-1]
            j -= 1
            
            
            
    # finish when meeting 0 on j or i
    if j != 0:
        while j > 0:
            
            aligned_sequence1 += seq1[j-1]
            aligned_sequence2 += '-'
            j -= 1
            
    
    aligned_sequence1 = aligned_sequence1[::-1]
    aligned_sequence2 = aligned_sequence2[::-1]
    
    return(aligned_sequence1, aligned_sequence2)

