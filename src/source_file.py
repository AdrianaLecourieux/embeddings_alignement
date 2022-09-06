# Source file


"""
with open("../data/fasta_sequences/5_3_EXONUCLEASE_1BGXT.fasta") as f:
    firstline = f.readline().rstrip()

print(firstline)
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

import numpy as np
"""from scipy.stats import pearsonr
corr= pearsonr(seq1, seq2)

print(corr)
"""
