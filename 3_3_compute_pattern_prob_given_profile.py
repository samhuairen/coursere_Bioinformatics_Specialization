profile = {
    'A': [0.2, 0.2, 0, 0, 0, 0, 0.9, 0.1, 0.1, 0.1, 0.3, 0],
    'C': [0.1, 0.6, 0, 0, 0, 0, 0, 0.4, 0.1, 0.2, 0.4, 0.6],
    'G': [0, 0, 1, 1, 0.9, 0.9, 0.1, 0, 0, 0, 0, 0],
    'T': [0.7, 0.2, 0, 0, 0.1, 0.1, 0, 0.5, 0.8, 0.7, 0.3, 0.4]
}


def compute_pattern_prob_by_given_profile(pattern, profile):
    prob = 1
    pat_len = len(pattern)
    for i in range(0, pat_len):
        nucleotide = pattern[i]
        prob_nucletide = profile[nucleotide][i]
        prob *= prob_nucletide
    return prob

out = compute_pattern_prob_by_given_profile('TCGTGGATTTCC', profile)