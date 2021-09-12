def immediate_neighbors(pattern):
    neighborhood = [pattern]
    for i in range(0, len(pattern)):
        nucleotide = pattern[i]
        for x in "ATCG".replace(nucleotide, ""):
            neighbor = pattern[:i] + x + pattern[i+1:]
            neighborhood.append(neighbor)
    return neighborhood

resu = immediate_neighbors("CAA")

resu2 = immediate_neighbors("CCAGTCAATG")