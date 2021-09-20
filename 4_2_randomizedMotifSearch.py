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


def produce_kmers(seq, k):
    kmers = []
    for i in range(0, len(seq)-k+1):
        window = seq[i:i+k]
        kmers.append(window)
    return kmers


def random_select_motifs(dna_list, k, t):
    from random import choices
    all_motifs = []
    for i in range(0, t):
        seq = dna_list[i]
        kmers = produce_kmers(seq, k)
        all_motifs.append(kmers)
    random_motifs = [choices(item, k=1) for item in all_motifs]
    return ["".join(item) for item in random_motifs]


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


def Motifs(profile, dna_list, k):
    motifs = []
    for i in range(0, len(dna_list)):
        motif = profile_most_probable_kmer(dna_list[i], k, profile)
        motifs.append(motif)
    return motifs
    

def RandomizedMotifSearch(seq_list, k, t):
    motifs_init = random_select_motifs(seq_list, k, t)
    best_motif = motifs_init
    while True:
        profile = profile_motifs(best_motif)
        motifs = Motifs(profile, seq_list, k)
        if score_motifs(motifs) < score_motifs(best_motif):
            best_motif = motifs
        else:
            break
    return best_motif


def select_best_motifs(dna_list, k, t):
    i = 0
    first_motifs = RandomizedMotifSearch(dna_list, k, t)
    while i < 1000:
        best_motif = RandomizedMotifSearch(dna_list, k, t)
        if score_motifs(first_motifs) > score_motifs(best_motif):
            first_motifs = best_motif[:]
        i += 1
    return first_motifs


seq_file = open("dataset_161_5-3.txt", "r")
k, t = map(int, seq_file.readline().strip().split(" "))
dna_list = []
sequences = seq_file.read()
for line in sequences.strip().split(" "):
    dna_list.append(line)
    
out3 = select_best_motifs(dna_list, k, t)
print(" ")
        
        
        
        
#################################
# my own way to find the minimum score motifs for all 1000 motifs results
best_motif_list = []
i = 0
while i < 1000:
    best = RandomizedMotifSearch(dna_list, k, t)
    best_motif_list.append(best)
    i += 1

scores_best = [score_motifs(motif) for motif in best_motif_list]
minimum = min(scores_best)
best_motifs = []
for i in range(0, len(scores_best)):
    if scores_best[i] == minimum:
        best_motifs.append(best_motif_list[i])

dna_list2 = ['TGACGTTC', 'TAAGAGTT', 'GGACGAAA', 'CTGTTCGC']
motif_init = ['TGA', 'GTT', 'GAA', 'TGT']
profile1 = profile_motifs(motif_init)
new = Motifs(profile1, dna_list2, 3)