def consensus(motifs):
    from collections import Counter
    consensus_motif = ''
    for i in range(0, len(motifs[0])):
        column_nucleotides = [motif[i] for motif in motifs]
        count_nucleotide = Counter(column_nucleotides)
        most_common = count_nucleotide.most_common()[0][0]
        consensus_motif += most_common
    return consensus_motif

