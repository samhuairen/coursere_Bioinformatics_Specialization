def hamming_distance(s1, s2):
    distance = 0
    for (i,j) in zip(s1, s2):
        if i != j:
            distance += 1
    return distance


def neighbor(pattern, d):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {"A", "T", "C", "G"}
    neighborhood = set()
    suffix_neighbor = neighbor(pattern[1:], d)
    for string in suffix_neighbor:
        if hamming_distance(string,pattern[1:]) < d:
            for nucleotide in "ATCG":
                neighborhood.add(nucleotide+string)
        else:
            neighborhood.add(pattern[0]+string)
    return neighborhood


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


def motifEnumeration(seq_set, k, d):
    patterns = set()
    kmers = []
    for seq in seq_set:
        for i in range(0, len(seq)-k+1):
            kmer = seq[i:i+k]
            kmers.append(kmer)
    for kmer in kmers:
        kmer_neighbor = neighbor(kmer, d)
        for pattern_ in kmer_neighbor:
            counts = [count_mismatch_pattern(seq, pattern_, d) for seq in seq_set]
            if all(count > 0 for count in counts):
                patterns.add(pattern_)
    return patterns


seq_set = ['CTTTTACGTAACAATGTGGCTCTAG', 'GGTCTGCTAGCCCTCTGTTCCTTGT', 'CGAAGTTCCCCTTAGTAGGCATAGA',
           'CGCTATCTAGACCGTCCTGAGACAT', 'CTGTGCATAGCGTGGCCGGAACCAA', 'AGTGTATTGAGCTAGGACGGAACCT']
k = 5
d = 2
out = motifEnumeration(seq_set,k,d)
