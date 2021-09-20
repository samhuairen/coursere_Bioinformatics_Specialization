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


def produce_weighted_kmer(seq, k, profile):
    from random import choices
    prob_kmer = []
    kmer_position = []
    for i in range(0, len(seq)-k+1):
        kmer = seq[i:i+k]
        prob = compute_pattern_prob_by_given_profile(kmer, profile)
        prob_kmer.append(prob)
        kmer_position.append(i)
    weight = [prob/sum(prob_kmer) for prob in prob_kmer]
    pos = choices(kmer_position, weights=weight, k=1)[0]
    return seq[pos:pos+k]


def gibbsSampler(dna_list, k, t, N):
    from random import randint
    motifs = random_select_motifs(dna_list, k, t)
    best_motif = motifs[:]
    for j in range(0, N):
        i = randint(0, t-1)
        motifs.pop(i)
        profile = profile_motifs(motifs)
        motif_i = produce_weighted_kmer(dna_list[i], k, profile)
        motifs.insert(i, motif_i)
        if score_motifs(motifs) < score_motifs(best_motif):
            best_motif = motifs
    return best_motif


def select_best_motifs(dna_list, k, t, N, S):
    i = 0
    first_motifs = gibbsSampler(dna_list, k, t, N)
    while i < S:
        best_motif = gibbsSampler(dna_list, k, t, N)
        if score_motifs(first_motifs) > score_motifs(best_motif):
            first_motifs = best_motif[:]
        i += 1
    return first_motifs

k = 8
t = 5
N = 100
S = 300
seq_list = ['CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
            'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
            'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']

out = select_best_motifs(seq_list, k, t, N, S)











seq_file = open("dataset_163_4-3.txt", "r")
k, t, N = map(int, seq_file.readline().strip().split(" "))
S = 20
dna_list = []
sequences = seq_file.read()
for line in sequences.strip().split(" "):
    dna_list.append(line)

out = select_best_motifs(dna_list, k, t, N, S)

###########
dna_list = []
seq_file = open("DosR.txt", "r")
sequences = seq_file.read()
for line in sequences.strip().split("\n"):
    dna_list.append(line)
    
t = 10
N = 1000
S = 100
for k in range(8, 20):
    out = select_best_motifs(dna_list, k, t, N, S)
    score = score_motifs(out)
    print(" ".join(out))
    print(score)
    

