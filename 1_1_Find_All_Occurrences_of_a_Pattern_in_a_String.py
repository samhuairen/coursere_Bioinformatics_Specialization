import re
file_download = open("rosalind_ba1d-3.txt", 'r')
pat = file_download.readline().strip()
sequence = file_download.read().strip()
find_results = re.finditer(pat, sequence)
out_pos = []

pat_len = len(pat)
seq_len = len(sequence)
for i in range(0, seq_len-pat_len +1):
    window = sequence[i:i+pat_len]
    if window == pat:
        out_pos.append(i)

out_file = open("out_results3.txt",'w')
for item in out_pos:
    out_file.write(str(item)+" ")
out_file.close()

##########


def find_all_occurences_of_a_string(seq, pattern):
    out_index = []
    len_seq = len(seq)
    len_pat = len(pattern)
    for i in range(0, len_seq-len_pat+1):
        window = seq[i:i+len_pat]
        if window == pattern:
            out_index.append(i)
    return out_index
print(find_all_occurences_of_a_string('GACGATATACGACGATA',"ATA"))