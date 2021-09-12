def hamming_distance(s1, s2):
    distance = 0
    for (i,j) in zip(s1, s2):
        if i != j:
            distance += 1
    return distance


def count_mismatch_pattern(seq, pattern, mismatch):
    count = 0
    len_seq = len(seq)
    len_pat = len(pattern)
    for i in range(0, len_seq - len_pat + 1 ):
        kmer = seq[i:i+len_pat]
        hm_distance = hamming_distance(pattern, kmer)
        if hm_distance <= mismatch:
            count += 1
    return count

seq = 'GGACTGGGAATGAATCAATTCGATTCGAACGTCCTTCCGTGCACGCCAGAATTCCACCGCGTGGGGTGCAGGACGGCTGCTCCGTGTTATTCTGCATCAACCAGGAAACTTCATGTAACCCACAGCTGGGAAGGCGAGTGCAATTCGCAGCTTGTGCAATGCCCCGTATGCTTTTTTTCGCGCATACCCACGCTCAGCAGGTAACGCAAAATTGCTAAAAGGGGTCATGGAGTCCCAATTACTCTTAACGAGTGGTCGCGACCCAGTCATATCGAGTCATACGTGATAGACACGGCGACCGGGGCAGCCGAGGATTCATCACTCGCCATCACAGGTTATTCGTTTTTTTCAGTCGTTAAAGGAGGGTCTTTGCATATATGCCAGGACTTGCAAT'
pat = "GTCGC"
mismatch = 3
out = count_mismatch_pattern(seq, pat, 3)
print(out)