def read_embedding(file):

    sequence = []
    with open(file, "r") as embedding:
        for line in embedding:
            vector = line.split()
            vector = [float(x) for x in vector]
            sequence.append(vector)
    return sequence

print(read_embedding("../5_3_exonuclease_1bgxt.t5emb"))