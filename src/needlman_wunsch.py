# --- Needlman and Wunsch alignment script

from pickle import FALSE # ??
import numpy as np

# --- With gap penalty fixed to 0

def transformation_NW(dot_matrix, seq1, seq2):
    
    # Creation of the transformed matrix
    
    col_seq1 = len(seq1) + 1
    row_seq2 = len(seq2) + 1
    transformed_matrix= np.zeros((col_seq1, row_seq2),dtype = int)
    for i in range(1,col_seq1):
        
            for j in range(1,row_seq2):
                
                top = transformed_matrix[i][j-1] 
                diagonal = transformed_matrix[i-1][j-1] + dot_matrix[i-1][j-1] 
                left = transformed_matrix[i-1][j] 
                transformed_matrix[i,j] = max(diagonal, top, left)
                
    return(np.matrix(transformed_matrix)) 

def needlman_wunsch(seq1,seq2, transformed_matrix):

    # Find optimal path
    
    aligned_sequence1 = ""
    aligned_sequence2 = ""
    i = len(seq1)
    j = len(seq2)

    #From the right bottom to the left top
    
    while i > 0 and j > 0:
        
        score_diagonal = transformed_matrix[i-1,j-1]
        score_top = transformed_matrix[i,j-1]
        score_left = transformed_matrix[i-1,j]
        max_value = max(score_diagonal, score_left, score_top)        

        # Calcule the score value
            
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





    
# --- With affine gap penalty set to -1 for gap opening and 0 for gap extension

def transformation_NW_affine_gap_penalty(dot_matrix, seq1, seq2):
    
    # Creation of the transformed matrix of zeros
    
    col_seq1 = len(seq1) + 1
    row_seq2 = len(seq2) + 1
    
    transformed_matrix= np.zeros((col_seq1, row_seq2),dtype = int)
   
   
    # Create a penalty matrix of 1 (no penalty) and gap penalty
    
    penalty_matrix = np.ones((col_seq1, row_seq2),dtype = int)
    open_gap = -1
    extension_gap = 0
    
    for i in range(0, col_seq1):
        
            for j in range(0, row_seq2):
                
                # First cell
                
                #print(i,j)
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
                    if max(top, left, diagonal) == left:
                        transformed_matrix[i,j] = left                        
                        if penalty_matrix[i, j-1] == 1:
                            penalty_matrix[i,j] = open_gap
                        
                        else: # sinon c'est extension de gap
                            penalty_matrix[i,j] = extension_gap
                    
                    elif max(top, left, diagonal) == top: # si le max c'est top et que 1 dans penlty alors ouverture
                        transformed_matrix[i,j] = top
                        if penalty_matrix[i-1, j] == 1:
                            penalty_matrix[i,j] = open_gap
                        
                        else: # sinon extension
                            penalty_matrix[i,j] = extension_gap
                    
                    else:
                        transformed_matrix[i][j] = diagonal
                        penalty_matrix[i,j] = 1                               
        # -------       
    return(np.matrix(transformed_matrix)) 