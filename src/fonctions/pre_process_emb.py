# --- Pre-processing of embeddings files

import numpy as np # import the numpy module and rename it
      
def read_embedding(file_emb):
    """Open and read embedding file.

    Split the lines to separate variables and transform them to float.
    Add a list per vector.

    Parameters
    ----------
    file : embedding file
        1024 variables per line for each vector
        
    Returns
    -------
    embedding_list : embedding list
        1024 variables per list into a major list.
    """
    embedding_list = []
    
    with open(file_emb, "r") as embedding:
        
        for line in embedding:
            
            vector = line.split()
            vector = [float(x) for x in vector]
            embedding_list.append(vector)  
          
    return(embedding_list)



 

    
