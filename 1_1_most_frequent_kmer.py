def find_most_frequent(seq, k):
    '''find the most frequent kmer with length k in a sequence '''
    from collections import defaultdict
    freq_table = defaultdict(int)
    seq_len = len(seq)
    for i in range(0, seq_len-k+1):
        window = seq[i:i+k]
        freq_table[window] += 1
    sorted_freqTable = sorted(freq_table.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
    return sorted_freqTable[0][0]


seq="TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT"
result=find_most_frequent(seq, k=3)