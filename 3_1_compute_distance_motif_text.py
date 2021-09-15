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


def find_motif_pattern_text(pattern, seq):
    """ find the kmer that minimize the distance between the patten and the kmer in seq"""
    len_pat = len(pattern)
    len_seq = len(seq)
    distance = []
    pattern_ = []
    for i in range(0, len_seq - len_pat +1):
        window = seq[i:i+len_pat]
        hm_dist = hamming_distance(window, pattern)
        distance.append(hm_dist)
        pattern_.append(window)
    distance_sorted = sorted(enumerate(distance), key=lambda x1: x1[1])
    min_distance_index = distance_sorted[0][0]
    return pattern_[min_distance_index]


def distanceBetweenPatternAndStrings(pattern, seq_set):
    k = len(pattern)
    distance = 0
    for seq in seq_set:
        dist = compute_distance_motif_text(pattern, seq)
        distance += dist
    return distance





seq_set = ['TTACCTTAAC', 'GATATCTGTC', 'ACGGCGTTCG',
           'CCCTAAAGAG', 'CGTCAGAGGT']
out = distanceBetweenPatternAndStrings("AAA", seq_set)
        

seq_file = open("dataset_5164_1-4.txt", 'r')
content = seq_file.readlines()
pattern = content[0].strip()
seq_set = [seq.strip() for seq in content[1].split(" ")]
for line in seq_file.readlines():
    seq_set.append(line.strip())
out = distanceBetweenPatternAndStrings(pattern, seq_set)