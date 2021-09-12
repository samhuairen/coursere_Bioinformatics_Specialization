def hamming_distance(s1, s2):
    distance = 0
    for (i,j) in zip(s1, s2):
        if i != j:
            distance += 1
    return distance


with open("dataset_9_3.txt", "r") as f:
    seqs = f.readlines()
    seq1 = seqs[0].strip()
    seq2 = seqs[1].strip()
    print(seq1)
    print(seq2)
    print(hamming_distance(seq1, seq2))