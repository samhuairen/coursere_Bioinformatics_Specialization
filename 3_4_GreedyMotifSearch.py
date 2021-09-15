motifs = [
    'ATTATATATAT',
    'ATGTGTATAGA',
    'ATAGATATAGA',
    'ATAGATAGATA'
]

motifs2 = [
    'TCGGGGGTTTTT',
    'CCGGTGACTTAC',
    'ACGGGGATTTTC',
    'TTGGGGACTTTT',
    'AAGGGGACTTCC',
    'TTGGGGACTTCC',
    'TCGGGGATTCAT',
    'TCGGGGATTCCT',
    'TAGGGGAACTAC',
    'TCGGGTATAACC'
]


def profile_motifs(motifs):
    out_profile = {"A": [], "C": [], "G": [], "T": []}
    for i in range(0, len(motifs[0])):
        column_nucleotides = [motif[i] for motif in motifs]
        a_frq = (column_nucleotides.count("A") + 1) / (len(column_nucleotides)+4)
        c_frq = (column_nucleotides.count("C") + 1) / (len(column_nucleotides)+4)
        g_frq = (column_nucleotides.count("G") + 1) / (len(column_nucleotides)+4)
        t_frq = (column_nucleotides.count("T") + 1) / (len(column_nucleotides)+4)
        out_profile['A'].append(a_frq)
        out_profile['C'].append(c_frq)
        out_profile['G'].append(g_frq)
        out_profile["T"].append(t_frq)
    return out_profile


def score_motifs(motifs):
    from collections import Counter
    scores = 0
    for i in range(0, len(motifs[0])):
        column_nucleotides = [motif[i] for motif in motifs]
        freq = Counter(column_nucleotides)
        unpopular_freq = len(column_nucleotides) - freq.most_common()[0][1]  # len of column minus the most popular freq
        scores += unpopular_freq
    return scores
    

def form_init_motif_matrix(seq_list, k):
    init_motif_list = []
    for i in range(len(seq_list)):
        seq = seq_list[i]
        init_motif_list.append(seq[0:k])
    return init_motif_list


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
        prob_nucleotide = profile[nucleotide][i]
        prob *= prob_nucleotide
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


def greedMotifsSearch(seq_list, k, t):
    """ """
    best_motifs = form_init_motif_matrix(seq_list, k)
    for i in range(0, len(seq_list[0])-k + 1):
        kmer = seq_list[0][i:i+k]
        motif1 = kmer
        motifs = [motif1]
        for j in range(1, t):
            profiles = profile_motifs(motifs)
            motif_i = profile_most_probable_kmer(seq_list[j], k, profiles)
            motifs.append(motif_i)
        if score_motifs(motifs) < score_motifs(best_motifs):
            best_motifs = motifs
    return best_motifs




seq_file = open("dataset_160_9.txt", "r")
k, t = map(int, seq_file.readline().strip().split(" "))
seq_list = []
sequences = seq_file.read()
for line in sequences.strip().split(" "):
    seq_list.append(line)
out = greedMotifsSearch(seq_list, k, t)

profile2 = {
    "A": [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
    "C": [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
    "G": [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
    "T": [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]
}

motif_list = ['AAGTGA','AGGTCA','ACGTTA','ACGTTT','ACGCGA','ATGCTA']

probs = [ compute_pattern_prob_by_given_profile(motif, profile2) for motif in motif_list]

seq_list_test = [
    'CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC',
    'GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC',
    'GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG',
]

