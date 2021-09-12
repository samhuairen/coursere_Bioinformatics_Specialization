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


def finding_clump(seq, k, L, t):
    patterns = set()
    seq_len = len(seq)
    for i in range(0, seq_len - L + 1):
        long_window = seq[i:i+L]
        fre_table = frequent_table(long_window, k)
        for item in fre_table.keys():
            if fre_table[item] >= t:
                patterns.add(item)
    return " ".join(patterns)


seq = "ATAGTCGCTAGGAGGAACGCTACAGGAACGCCACGGGCTTAAAATCGCCTCGATCAAATTGAGGATGGCGATACAAATTACACAAATTATAATTAGAGGCCACGTGTGCCGATAAGATACTCAGCTTCCGACAGCACCAGAGTTCAACTTCAAGCGAATAATTTGGTCAGCCTGATGTTAAAATGATGGGCCCTTCTACACGCCTATGGTGTCTGGTTTCAACTGATGCTTCCTGCAGGTGGGAGACTACGGGAGCTCCCACGTGCCTTCCCCTGGGAAGTGGAGGCGCTGCGATGTTGACCGAGTTAGGCCGATTGAGCGTAGCGGTATTAGTCATTAGTGTGTAGTATTAGTCTGCTACTTTACCTGGCGCGAATGCACACCCGTGCACTAGAAAGGTTACACGACTAGCCGAGATGACGAATTGCGGTATAATAGCCGACCCCTGCGGCCCATCTCTTAGTACTCTGCCCACTGTACCTCCATCACAGCAAAAGACTTGCGACCTGTAGTCATCCTTGATGGGTAAGAGGGATACTCATGGTGGTATCAGTGGCCGCGCGAACACCTTACAAGCCCACTATACGCTTAACCAGATTGGCGTTTACCCGCCTGTGCCGTGCCGTGGGGTGGGGCCGTGGGGAACGTAGGTTGGGGGCCGAGAGCTGCTCTAACTCTCACATAACGACTCTGTGACGCTCCGCCTGTCTGTGCATGTATTCAAGTCTACACACAATGGAAAAGAATGAGGCCACTACCCAGCGATATATGAAGCGAGGCATGTTAACCCGAGTTTTGCGTGGACTCAAGAGCAAGAGGACTCAAGTAAAAAAAAAAAAGGAAAAAAGGGGGGTATTAATGGCCTTCAGAGCGTAGGTCCTATCGAATCGAAATCGTCGAAAAAAAACCGCAGTCACCGCCGCAAGCAGAGTCCAGAGTCCTTAGAAGAGTCCTCTTTTTTACAGTTATGATTAGGGGCGACTGCTATCGGTCAGCACTCGAGGGTTTCCACCTAGTTGTCAGGCGATCTCGCTGGAGGACGGTATGAACGAGTACGTATAATCGCCCCTCTCCCACCGGCATTCAGCCGTGGTCCGCTGAAATATGTCCAGTCCGGGAAAGAACGCCTCGTTACGCTGATACAGAGATAAATCCGGGCGGGTGCTCACCGTTTAGGAAACAAACCCACTGAAGGAACGACCGACTTATGTATGAAAGGGAACCATACCATGCAGGCTAGTGAAGGGTATGAGCTTCCCCTAGTGGTTGTTGGGAGAGCGTAGAGCTGGACCCAAGGTTATCAACTCCATCTCGCAACTCGGCCCTCAGGCCCCCAGCGCGCAGTCTCTAGCATCCGGTTCACTACTTACCGGAAGTAACCTATGACTGGGTTCGCTAACCCCAGAAAGAACATTACCAGAGGGGATGGGATGGTTGGATGGTTGCAAACAGAATTACGCGAATAGAAACATGGTGTGGAGTGTGGAGGTGTGGAGATAAGGGCCTTAATGATCCAAATAATCAATAGCTAGTAAGTTCGAGTACTTTAGACAGGGCCTACAAATTATGTGCCGTTCGCGTAAATCTAGAATTAGTATAGTAGGTCTTCAATTCTCGTAATAGAAAGCTGAATTCCGGCCCCCTACAAGTCAATACTCGCATTGGATACTAACGATAGGCACCGTACGGGGAGCTCCTTATATTCAAACGAAAGGAAAGAACCAGAAACGAAAGAGGAAAAGACCGTAAAATGCGGCGCTCATCATGACTGCTTAACTTGGCTATTGTGCGCGAGAACGTGGCTTTGGGGTCTATCGAGGTGGAGTACCATAGTGCGACTCGTTTAGCACTCCAGCACTGGAGTACTTTTAACTCTGGACACAGAATCCTTGCCTAGTAGTTCAAGTTGAGCGGAGACGTTTCGCTGGACGTTTCTTTCTCTTCTCCTTCGGTGC"
k = 8
L = 30
t = 3

result = finding_clump(seq, k, L, t)

seq = open("/Users/zhanghuairen/PycharmProjects/unite_python3/E_coli.txt").read().strip()
k = 9
L = 500
t = 3
result = finding_clump(seq, k, L, t)