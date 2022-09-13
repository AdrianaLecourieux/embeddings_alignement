# --- Creation of the score matrix with calculating dot product

import numpy as np # import the numpy module and rename it


def dot_product(emb1_list, emb2_list):
    """Open and read embedding file.

    Calculation of the dot product obtained 
    from the embedding vectors of one position of protein 1 and another
    position of protein 2. 
    A transposition of one list is necessary to calculate the dot product

    Parameters
    ----------
    emb1_list : embeddings list
        1024 variables per line for each vector
    emb2_list : embeddings list
        1024 variables per line for each vector
        
    Returns
    -------
    calcul_dot_product : numpy array
        An array where each value corresponds to the dot product obtained 
        from the embedding vectors of one position of protein 1 and another
        position of protein 2.
    """
    calcul_dot_product = np.dot(emb2_list, np.array(emb1_list).T) 

    return(calcul_dot_product)
  