def reverse_complement(seq):
    base_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    seq_complement = ""
    for i in seq:
        seq_complement += base_dict[i]
    return seq_complement[::-1]


def hamming_distance(s1, s2):
    distance = 0
    for (i, j) in zip(s1, s2):
        if i != j:
            distance += 1
    return distance


def count_mismatch_pattern(seq, pattern, mismatch):
    count = 0
    len_seq = len(seq)
    len_pat = len(pattern)
    for i in range(0, len_seq - len_pat + 1):
        kmer = seq[i:i + len_pat]
        hm_distance = hamming_distance(pattern, kmer)
        if hm_distance <= mismatch:
            count += 1
    return count


def find_frequent_words_with_mismatches_and_reversed(seq, k, d):
    """ find  all k-mers Pattern maximizing the sum Countd(Text, Pattern)+ Countd(Text, Patternrc)
     over all possible k-mers"""
    from itertools import product
    out_pairs = []
    freq_dict = {}
    nucleotides = ["A", "T", "C", "G"]
    patterns_forward = ["".join(item) for item in product(nucleotides, repeat=k)]
    patterns_reverse = [reverse_complement(pat) for pat in patterns_forward]
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


with open('Salmonella_enterica.txt', 'r') as f:
    seq = ""
    for line in f:
        if line.startswith(">"):
            continue
        else:
            seq += line.strip()

skew_position = find_minimum_skew(seq)
seq_analysis1 = seq[3764856-250: 3764858+250]
seq_analysis2 = seq[3764856-500: 3764858+500]

out1 = find_frequent_words_with_mismatches_and_reversed(seq_analysis1, 9, 1)
out2 = find_frequent_words_with_mismatches_and_reversed(seq_analysis2, 9, 1)

s1 = "TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC"
s2 = 'GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA'

hm = hamming_distance(s1, s2)
print(hm)

s3 = 'GATACACTTCCCGAGTAGGTACTG'
s3_out = find_minimum_skew(s3)

s4 ='TACGCATTACAAAGCACA'
pat = 'AA'
s4_out = count_mismatch_pattern(s4, pattern=pat, mismatch=1)

