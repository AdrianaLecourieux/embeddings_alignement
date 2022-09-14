# --- Needleman and Wunsch alignment script

from pickle import FALSE, TRUE
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
    seq1: list
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
  
    
    # we want the begin with the best score in the last column 
    # the max of i in j col
    col = len(seq1)
    index_max = np.where(transformed_matrix[:,col] == np.amax(transformed_matrix[:,col]))
    
    
    idx_row = index_max[0]
    
   
    align_list = []
    
    for i in idx_row: 
        
        j = col
          # Find optimal path
    
        aligned_sequence1 = ""
        aligned_sequence2 = ""
        
        check_zero = FALSE # check index
        
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
                
            if max_value == 0:
                check_zero = TRUE
                break
        
    # finish when index in j = 0
    
        if not check_zero:
    
            if j != 0:
                while j > 0:
                    
                    aligned_sequence1 += seq1[j-1]
                    aligned_sequence2 += '-'
                    j -= 1
                    
    
        aligned_sequence1 = aligned_sequence1[::-1]
        aligned_sequence2 = aligned_sequence2[::-1]
        
        align_list.append([aligned_sequence1, aligned_sequence2])
        
        
    return(align_list)

def transformation_semi_global_gp(dot_matrix, seq1, seq2):
    """Creation of the transformed matrix from the dot product matrix and 
    fasta sequences.

    The transformed matrix is created by following the semi-global algorithm (based
    on NW algorithm) 
    with the dot product matrix as the score matrix and a affine gap penalty set to:
    -1 for gap opening and 0 for gap extension

    Parameters
    ----------
    dot_matrix : array
        An array where each value corresponds to the dot product obtained 
        from the embedding vectors of one position of protein 1 and another
        position of protein 2.
    seq1: list
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
  
    # Create a penalty matrix of 1 (no penalty) and gap penalty
    penalty_matrix = np.ones((row_seq, col_seq),dtype = int)
    open_gap = -1
    extension_gap = 0
    
    for i in range(0,row_seq):
        
            for j in range(0,col_seq):
                # First cell
                
                if i == 0 and j == 0:
                    transformed_matrix[i][j] = 0
                
                # First row
                
                elif i == 0:
                     
                    if penalty_matrix[i, j-1] == 1 : # no penalty
                        
                        transformed_matrix[i][j] = transformed_matrix[i][j-1] 
                        
                    else:
                        
                        transformed_matrix[i][j] = transformed_matrix[i][j-1] + penalty_matrix[i, j-1] # penalty
                        
                    # Update penalty matrix
                     
                    if penalty_matrix[i, j-1] == 1: 
                        penalty_matrix[i,j] = open_gap
                    
                    else:
                        penalty_matrix[i,j] = extension_gap
                
                # First column
                
                elif j == 0:
                   
                    if penalty_matrix[i-1, j] == 1:
                        transformed_matrix[i][j] = transformed_matrix[i-1][j]
                    else:
                        transformed_matrix[i][j] = transformed_matrix[i-1][j] + penalty_matrix[i-1, j]
                    
                     # Update penalty matrix
                     
                    if penalty_matrix[i-1, j] == 1:
                        penalty_matrix[i,j] = open_gap
                    
                    else:
                        penalty_matrix[i,j] = extension_gap
                
                # Rest of the matrix
                
                else:
        # ------- Pour avoir le score autour de i,j
                # si pas de penalité dans case gauche alors juste le score de gauche
                    if penalty_matrix[i, j-1] == 1 :
                         left= transformed_matrix[i][j-1] 
                    else: # sinon score + penalité
                        left = transformed_matrix[i][j-1] +  penalty_matrix[i, j-1]
                        
                    # si pas de penalité au top ..... score du haut
                    if penalty_matrix [i-1][ j] == 1 :
                        
                        top = transformed_matrix[i-1][j] 
                    else: # sinon ... + penalité
                        top =transformed_matrix[i-1][j] + penalty_matrix [i-1][ j]
                        
                       # pas de penalité 
                    diagonal = transformed_matrix[i-1][j-1] + dot_matrix[i-1][j-1]
        # -------
        # ------- Check if it's gap opening or gap extension          
                    #si max c'est left alors prend left et si c'est 1 dans la penality matrix alors c'est ouverture de gap
                    max_val = max(top, left, diagonal)
                    if max_val == left:
                        transformed_matrix[i,j] = left                        
                        if penalty_matrix[i, j-1] == 1:
                            penalty_matrix[i,j] = open_gap
                        
                        else: # sinon c'est extension de gap
                            penalty_matrix[i,j] = extension_gap
                    
                    elif max_val == top: # si le max c'est top et que 1 dans penlty alors ouverture
                        transformed_matrix[i,j] = top
                        if penalty_matrix[i-1, j] == 1:
                            penalty_matrix[i,j] = open_gap
                        
                        else: # sinon extension
                            penalty_matrix[i,j] = extension_gap
                    
                    else:
                        transformed_matrix[i][j] = diagonal
                        penalty_matrix[i,j] = 1             

            
    return(np.matrix(transformed_matrix)) 