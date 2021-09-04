
def count_string(seq, string):
    times = 0
    len_seq = len(seq)
    len_string = len(string)
    for i in range(0, len_seq-len_string+1):
        window = seq[i:i+len_string]
        if window == string:
            times += 1
    return times


def frequent_table(seq, k):
    len_seq = len(seq)
    fre_table = {}
    for i in range(0, len_seq-k+1):
        window = seq[i:i+k]
        if window not in fre_table.keys():
            fre_table[window] = 1
        else:
            fre_table[window] += 1
    return fre_table


def frequent_patterns(seq, k):
    patterns = []
    fre_table = frequent_table(seq, k)
    print(fre_table)
    max_value = max(fre_table.values())
    for item in fre_table.keys():
        if fre_table[item] == max_value:
            patterns.append(item)
    return patterns


seq = 'ATCATCCAGATGGGGGGGCGCTGACGTCTCATCGATCATCCAGTCTCATCGTCCCTTATCGTCTCATCGATCATCCAGATGGGGGGTCCCTTATCGTCTCATCGGATGGGGGGGTCTCATCGGTCTCATCGGTCTCATCGGATGGGGGGGATGGGGGGATCATCCAGCGCTGACGATGGGGGGGATGGGGGGATCATCCAGTCTCATCGGTCTCATCGTCCCTTATCGATGGGGGGATCATCCATCCCTTATCATCATCCAGATGGGGGGATCATCCATCCCTTATCGTCTCATCGGATGGGGGGGATGGGGGGATCATCCATCCCTTATCGCGCTGACGCGCTGACGATGGGGGGGATGGGGGGATCATCCAGCGCTGACGCGCTGACGTCTCATCGTCCCTTATCGTCTCATCGATCATCCAGATGGGGGGGATGGGGGGGATGGGGGGATCATCCAGCGCTGACATCATCCATCCCTTATCGTCTCATCGGTCTCATCGGCGCTGACGATGGGGGGGTCTCATCGGTCTCATCGGATGGGGGGGTCTCATCGGTCTCATCGTCCCTTATCATCATCCAGTCTCATCGTCCCTTATCGCGCTGACTCCCTTATCGCGCTGACTCCCTTATCGTCTCATCGTCCCTTATCGTCTCATCGGCGCTGACGATGGGGGGTCCCTTATCATCATCCAATCATCCAGCGCTGACTCCCTTATCATCATCCAATCATCCAGTCTCATCGATCATCCAATCATCCAGTCTCATCGGTCTCATCGGCGCTGACGTCTCATCGGTCTCATCGGCGCTGACGATGGGGGGGTCTCATCGGATGGGGGG'
k = 11

results = frequent_patterns(seq, k)
print(results)