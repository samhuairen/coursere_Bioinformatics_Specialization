profile = {
    'A': [0.2, 0.2, 0.3, 0.2, 0.3],
    'C': [0.4, 0.3, 0.1, 0.5, 0.1],
    'G': [0.3, 0.3, 0.5, 0.2, 0.4],
    'T': [0.1, 0.2, 0.1, 0.1, 0.2]
}


def compute_pattern_prob_by_given_profile(pattern, profile):
    prob = 1
    pat_len = len(pattern)
    for i in range(0, pat_len):
        nucleotide = pattern[i]
        prob_nucletide = profile[nucleotide][i]
        prob *= prob_nucletide
    return prob


def profile_most_probable_kmer(seq, k, profile):
    prob_list = []
    len_seq = len(seq)
    for i in range(0, len_seq - k + 1):
        kmer = seq[i:i + k]
        prob = compute_pattern_prob_by_given_profile(kmer, profile)
        prob_list.append((kmer, prob))
    
    maximum_prob = max([item[1] for item in prob_list])
    for item in prob_list:
        if item[1] == maximum_prob:
            most_kmer = item[0]
            break
    return most_kmer


data_file = open("dataset_159_3-2.txt", "r")
seq = data_file.readline().strip()
k = int(data_file.readline().strip())
profile_dict = {}
profile_dict['A'] = [float(item) for item in data_file.readline().strip().split(" ")]
profile_dict['C'] = [float(item) for item in data_file.readline().strip().split(" ")]
profile_dict['G'] = [float(item) for item in data_file.readline().strip().split(" ")]
profile_dict['T'] = [float(item) for item in data_file.readline().strip().split(" ")]
data_file.close()

out = profile_most_probable_kmer(seq, k, profile_dict)