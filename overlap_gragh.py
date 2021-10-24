from collections import defaultdict


def overlap_gragh(string_list):
    adjacency_list = defaultdict(set)
    for pattern in string_list:
        adjacency_list[pattern[:-1]].add(pattern)

    for pattern in string_list:
        suffixes = adjacency_list[pattern[1:]]
        if suffixes:
            print(pattern+' -> ' + ",".join(adjacency_list[pattern[1:]]))
    return adjacency_list


data = open("dataset_198_10.txt", 'r')
string_list = data.read().strip().split("\n")
data.close()

overlap_gragh(string_list)


