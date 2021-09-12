def hamming_distance(s1, s2):
    distance = 0
    for (i,j) in zip(s1, s2):
        if i != j:
            distance += 1
    return distance


def neighbor(pattern, d):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {"A", "T", "C", "G"}
    neighborhood = set()
    suffix_neighbor = neighbor(pattern[1:], d)
    for string in suffix_neighbor:
        if hamming_distance(string,pattern[1:]) < d:
            for nucleotide in "ATCG":
                neighborhood.add(nucleotide+string)
        else:
            neighborhood.add(pattern[0]+string)
    return neighborhood