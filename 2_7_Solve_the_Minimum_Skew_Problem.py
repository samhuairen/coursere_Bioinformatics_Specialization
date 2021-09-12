

def find_minimum_skew(seq):
    out_positions = []
    len_seq = len(seq)
    skew_genome = []
    skew0 = 0
    for i in range(0, len_seq):
        if seq[i] == "G":
            skew_new = skew0 + 1
            skew_genome.append(skew_new)
        elif seq[i] == "C":
            skew_new = skew0 - 1
            skew_genome.append(skew_new)
        else:
            skew_new = skew0
            skew_genome.append(skew_new)
        skew0 = skew_new
    skew_genome.insert(0, skew0)
    minimum = min(skew_genome)
    for i in range(0, len(skew_genome)):
        if skew_genome[i] == minimum:
            out_positions.append(i)
    return out_positions

