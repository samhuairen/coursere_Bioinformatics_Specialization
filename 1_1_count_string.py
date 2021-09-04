def count_string(seq, string):
    times = 0
    len_seq = len(seq)
    len_string = len(string)
    for i in range(0, len_seq-len_string+1):
        window = seq[i:i+len_string]
        if window == string:
            times += 1
    return times


count_string('GACCATCAAAACTGATAAACTACTTAAAAATCAGT', 'AAA')