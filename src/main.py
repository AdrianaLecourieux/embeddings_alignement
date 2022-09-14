# --- Main File --- #


# --- Import modules

import argparse  # to interact with arguments

# --- Import functions from fonctions folder
from fonctions import pre_process_emb 
from fonctions import score_fonction
from fonctions import fasta_sequences 
from fonctions import needleman_wunsch as NW
from fonctions import smith_waterman as SW
from fonctions import semi_global_alignment as SG
from fonctions import save_output

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    
    # --- Required arguments
    
    parser.add_argument("-emb1","--embedding1",  help= "Enter Embedding 1 in .t5emb extension")
    parser.add_argument("-emb2","--embedding2",  help= "Enter Embedding 2 in .t5emb extension")
    parser.add_argument("-f1","--fasta1",  help= "Enter Fasta 1 in .FASTA extension")
    parser.add_argument("-f2","--fasta2",  help= "Enter Fasta 2 in .FASTA extension")
    
    # --- Optionnal argument
    
        # Method
    parser.add_argument("-m","--method",  help='Choose a "global" (Needleman and Wunsch), "local" (Smith and Waterman) or "semi_global" alignment algorithm.  -m global default"')
    
        # Gap penalty
    parser.add_argument("-g","--gap_penalty",  help='Use this option to add affine gap penalty (Enter "affine" to used -1 for gap opening and 0 for gap extension). Else, gap penalty is fixed to 0')
    
    # --- Assign arguments to variables
    
    args = parser.parse_args()    
    seq1 = args.embedding1
    seq2 = args.embedding2
    fasta1 = args.fasta1
    fasta2 = args.fasta2
    alignment_method = args.method
    gap_penalty = args.gap_penalty

    # Check the file extension
     
    if not seq1.endswith(".t5emb") or not seq2.endswith(".t5emb"):   
        raise Exception("ERROR : The two first arguments must be .t5emb extension")
    
    
    else: 
        embedding1_list = pre_process_emb.read_embedding(seq1)
        embedding2_list = pre_process_emb.read_embedding(seq2) 
        fasta1_list = fasta_sequences.read_fasta(fasta1)
        fasta2_list = fasta_sequences.read_fasta(fasta2)
        
                
        # Check the correspondance between embedding and fasta file 
        
        if len(embedding1_list) != len(fasta1_list) or len(embedding2_list) != len(fasta2_list):
            raise Exception("ERROR : For one sequence, embedding and fasta sequence must have the same length.")

        # Calcul of the dot product matrix
        
        dot_pro_mat = score_fonction.dot_product(embedding1_list, embedding2_list)
            
            # ------ Realise needleman and Wunsh alignment 
        
        if alignment_method == "global" or not alignment_method :
            
            # With affine gap penalty
            
            if gap_penalty == "affine":
                
                transformed_matrix_gp = NW.transformation_NW_affine_gap_penalty(dot_pro_mat, fasta1_list, fasta2_list)
                
                seq_aligned_list = NW.needleman_wunsch(fasta1_list,fasta2_list, transformed_matrix_gp)
                
                # Save output in .txt file
                
                output = save_output.save_in_txt_gp(seq1, seq2, alignment_method, seq_aligned_list)
                    
            # With fixed gap penalty
            
            else: 
                
                transformed_matrix = NW.transformation_NW(dot_pro_mat, fasta1_list, fasta2_list)
                
                seq_aligned_list = NW.needleman_wunsch(fasta1_list,fasta2_list, transformed_matrix)
                       
                # Save output in .txt file
                
                output = save_output.save_in_txt(seq1, seq2, alignment_method, seq_aligned_list)
        
        
             # ------ Realise Smith and Waterman alignment 
            
        elif alignment_method == "local":
            
            if gap_penalty == "affine":
                
                transformed_matrix_gp = SW.transformation_SW_gp(dot_pro_mat, fasta1_list, fasta2_list)
                seq_aligned_list = SW.smith_waterman(fasta1_list,fasta2_list, transformed_matrix_gp)
                
                # Save output in .txt file
                
                output = save_output.save_in_txt_gp(seq1, seq2, alignment_method, seq_aligned_list)
            
            else:
                
                transformed_matrix_SW = SW.transformation_SW(dot_pro_mat, fasta1_list, fasta2_list)
                seq_aligned_list = SW.smith_waterman(fasta1_list,fasta2_list, transformed_matrix_SW)
                
                # Save output in .txt file
                
                output = save_output.save_in_txt(seq1, seq2, alignment_method, seq_aligned_list)
            


            # ------ Realise Semi-global alignment 
        
        elif alignment_method == "semi_global":
            
            if gap_penalty == "affine":
                transformed_matrix_SW = SG.transformation_semi_global_gp(dot_pro_mat, fasta1_list, fasta2_list)
                seq_aligned_list = SG.semi_global(fasta1_list,fasta2_list, transformed_matrix_SW)   
                
                # Save output in .txt file
                  
                output = save_output.save_in_txt_gp(seq1, seq2, alignment_method, seq_aligned_list)
            else:
                transformed_matrix_SW = SG.transformation_semi_global(dot_pro_mat, fasta1_list, fasta2_list)
                seq_aligned_list = SG.semi_global(fasta1_list,fasta2_list, transformed_matrix_SW)     
                
                # Save output in .txt file
                
                output = save_output.save_in_txt(seq1, seq2, alignment_method, seq_aligned_list)