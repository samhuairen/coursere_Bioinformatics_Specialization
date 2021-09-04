
def show_all_start_pos_of_pattern(seq, pattern):
    index = []
    len_seq = len(seq)
    len_pat = len(pattern)
    for i in range(0, len_seq-len_pat+1):
        window = seq[i:i+len_pat]
        if window == pattern:
            index.append(i)
    return index

seq = open("/Users/zhanghuairen/PycharmProjects/unite_python3/Vibrio_cholerae.txt").read().strip()
pattern = 'CTTGATCAT'
index = show_all_start_pos_of_pattern(seq, pattern)