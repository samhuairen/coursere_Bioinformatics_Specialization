def hamming_distance(s1, s2):
    distance = 0
    for (i,j) in zip(s1, s2):
        if i != j:
            distance += 1
    return distance


def compute_distance_motif_text(pattern, seq):
    distance = []
    len_pat = len(pattern)
    len_seq = len(seq)
    for i in range(0, len_seq - len_pat + 1):
        window = seq[i:i+len_pat]
        hm_dist = hamming_distance(pattern, window)
        distance.append(hm_dist)
    return min(distance)


def distanceBetweenPatternAndStrings(pattern, seq_set):
    distance = 0
    for seq in seq_set:
        dist = compute_distance_motif_text(pattern, seq)
        distance += dist
    return distance


def medianString(seq_set, k):
    import math
    from itertools import product
    distance = math.inf
    all_patterns = ["".join(item) for item in product("ATCG", repeat=k)]
    for i in range(0, len(all_patterns)):
        pattern = all_patterns[i]
        if distance > distanceBetweenPatternAndStrings(pattern, seq_set):
            distance = distanceBetweenPatternAndStrings(pattern, seq_set)
            median = pattern
    return median

seq_set = ['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG', 'GCTGAGCACCGG', 'AGTTCGGGACAG']
k = 3
out = medianString(seq_set,k)

seq_file = open("dataset_158_9.txt", 'r')
content = seq_file.readlines()
k = int(content[0].strip())
seq_set = [seq.strip() for seq in content[1].split(" ")]
for line in seq_file.readlines():
    seq_set.append(line.strip())
out = medianString(seq_set, k)