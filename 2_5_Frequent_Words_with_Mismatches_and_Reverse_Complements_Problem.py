
def reverse_seq(seq):
    base_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    seq_reverse = ""
    for i in seq:
        seq_reverse += base_dict[i]
    return seq_reverse[::-1]


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


def find_frequent_words_with_mismatches(seq, k, d):
    """ find the frequent patterns length k with mismatches less than d nucleotides"""
    from itertools import product
    out_patterns = []
    freq_dict = {}
    nucleotides = ["A", "T", "C", "G"]
    all_patterns = ["".join(item) for item in product(nucleotides, repeat=k)]
    for pat in all_patterns:
        count = count_mismatch_pattern(seq, pat, d)
        if count not in freq_dict.keys():
            freq_dict[count] = [pat]
        else:
            freq_dict[count].append(pat)
    max_pattern_number = max(freq_dict)
    for key, value in freq_dict.items():
        if key == max_pattern_number:
            out_patterns.append(value)
    return out_patterns


def find_frequent_words_with_mismatches_and_reversed(seq, k, d):
    """ find  all k-mers Pattern maximizing the sum Countd(Text, Pattern)+ Countd(Text, Patternrc)
     over all possible k-mers"""
    from itertools import product
    out_pairs = []
    freq_dict = {}
    nucleotides = ["A", "T", "C", "G"]
    patterns_forward = ["".join(item) for item in product(nucleotides, repeat=k)]
    patterns_reverse = [reverse_seq(pat) for pat in patterns_forward]
    pattern_pairs = [(pat, reverse_pat) for pat, reverse_pat in zip(patterns_forward, patterns_reverse)]
    
    for pat in pattern_pairs:
        count1 = count_mismatch_pattern(seq, pat[0], d)
        count2 = count_mismatch_pattern(seq, pat[1], d)
        count = count1 + count2
        if count not in freq_dict.keys():
            freq_dict[count] = [pat]
        else:
            freq_dict[count].append(pat)
    max_pattern_number = max(freq_dict)
    for key, value in freq_dict.items():
        if key == max_pattern_number:
            out_pairs.append(value)
    return out_pairs


seq = 'CCCCGTGCCCAGGTGCCGTGCCCCCAGCCCCCCCCCCCCCCGTGCAGGACAGCCGTGGAGACAGCAGGTGCCGACCCCCCCAGCCCAGCCGACCCCCAGCCCCCCGACCCCCCCCCAGCAGCCCCCAGCCCAGGACCCCCAGGTGCCCCCCCCCCCCCCGTGCAGCCCCCCGTGCCCCGTGCCCAGGACCGACAGCAGCCCAGCCGTGCAGGTG'
k = 6
d = 2
out = find_frequent_words_with_mismatches_and_reversed(seq, k, d)