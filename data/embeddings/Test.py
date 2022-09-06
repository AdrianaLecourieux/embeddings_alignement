"""import argparse

parser = argparse.ArgumentParser()
parser.parse_args()
"""

"""
embeddings = open("5_3_exonuclease_1bgxt.t5emb", "r")

print(embeddings.readline(1024))

fd = open('5_3_exonuclease_1bgxt.t5emb', 'r')
n = 0
while fd.readline():
    n += 1
    print(n)
"""

with open("5_3_exonuclease_1bgxt.t5emb") as f:
    firstline = f.readline().rstrip()

print(firstline)

#Une ligne est bien égale à un vecteur 