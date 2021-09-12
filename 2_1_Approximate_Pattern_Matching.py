def hamming_distance(s1, s2):
    distance = 0
    for (i,j) in zip(s1, s2):
        if i != j:
            distance += 1
    return distance


def find_approximate_pattern_pos(seq, pattern, mismatch_number):
    out_position = []
    len_seq = len(seq)
    pat_len = len(pattern)
    for i in range(0, len_seq - pat_len+1):
        kmer = seq[i:i+pat_len]
        hm_distance = hamming_distance(pattern, kmer)
        if hm_distance <= mismatch_number:
            out_position.append(i)
    return out_position


with open("dataset_9_4-2.txt", 'r') as f:
    pattern = f.readline().strip()
    seq = f.readline().strip()
    mismatch_number = int(f.readline().strip())
    resu = find_approximate_pattern_pos(seq, pattern, mismatch_number)
    print(" ".join([ str(item) for item in resu]))