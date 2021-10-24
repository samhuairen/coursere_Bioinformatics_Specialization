from collections import defaultdict
from itertools import product

binary_strings = ["".join([str(number) for number in item]) for item in list(product((0, 1), repeat=4))]


def overlap_gragh(string_list):
    adjacency_list = defaultdict(set)
    for pattern in string_list:
        adjacency_list[pattern[:-1]].add(pattern)

    for pattern in string_list:
        suffixes = adjacency_list[pattern[1:]]
        if suffixes:
            print(pattern+' -> ' + ",".join(adjacency_list[pattern[1:]]))
    return


