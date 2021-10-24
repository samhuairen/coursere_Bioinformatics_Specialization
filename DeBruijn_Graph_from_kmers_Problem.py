from collections import defaultdict


def DeBruijn(patterns):
    adjecency_list = {}
    for pattern in patterns:
        if pattern[:-1] not in adjecency_list:
            adjecency_list[pattern[:-1]] = []
        adjecency_list[pattern[:-1]].append(pattern[1:])
    
    adjecency_dict_sort = dict(sorted(adjecency_list.items()))
    return adjecency_dict_sort

def DeBruijn(patterns):
    adjecency_list = defaultdict(list)
    for pattern in patterns:
        adjecency_list[pattern[:-1]].append(pattern[1:])
    return adjecency_list

#kmers = ['GAGG',  'CAGG', 'GGGG', 'GGGA', 'CAGG', 'AGGG', 'GGAG']
#graphs = DeBruijn(kmers)
#for key, value in graphs.items():
#print(key + ' -> ' + ",".join(value) )




data = open('dataset_200_8-6.txt', 'r')
patterns = data.read().strip().split("\n")
data.close()

with open('out_adjencency_list6.txt', 'w') as f:
    graphs = DeBruijn(patterns)
    for key, value in graphs.items():
        f.write(key + ' -> '+",".join(value) + "\n")

