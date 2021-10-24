
def generate_edges(text, k):
    edges = []
    for i in range(0, len(text)-k+1):
        kmer = text[i:i+k]
        edges.append(kmer)
    return edges


def generate_nodes(text, k):
    """ node is the overlap between consecutive edges, and for de Bruijn graph, nodes
    are not repetive"""
    nodes = set()
    for i in range(0, len(text)-(k-1)+1):
        kmer = text[i:i+k-1]
        nodes.add(kmer)
    return nodes


def deBruijn_graph(text, k):
    adjecency_list = {}
    edges = generate_edges(text, k)
    nodes = generate_nodes(text, k)
    for node in nodes:
        suffes = []
        for edge in edges:
            if node == edge[:-1]:
                suffes.append(edge[1:])
        adjecency_list[node] = suffes
    return adjecency_list


data = open('dataset_199_6-5.txt', 'r')
k = int(data.readline().strip())
text = data.readline().strip()
data.close()

adjecency_list = deBruijn_graph(text, k)
with open("adjecency_list6.txt", 'w') as f:
    for key, value in adjecency_list.items():
        if value:
            f.write(key + ' -> ' + ','.join(value)+"\n")


adjecency_list = deBruijn_graph('TAATGCCATGGGATGTT', 2)
for key, value in adjecency_list.items():
    if value:
        print(key + ' -> ' + ','.join(value))
        
adjecency_list2 = deBruijn_graph('TAATGCCATGGGATGTT', 3)
for key, value in adjecency_list2.items():
    if value:
        print(key + ' -> ' + ','.join(value))

adjecency_list3 = deBruijn_graph('TAATGCCATGGGATGTT', 4)
for key, value in adjecency_list3.items():
    if value:
        print(key + ' -> ' + ','.join(value))
        
text1 = 'TAATGCCATGGGATGTT'
text2 = 'TAATGGGATGCCATGTT'

adjecency_list_text1 = deBruijn_graph(text1, 3)
for key, value in adjecency_list_text1.items():
    if value:
        print(key + ' -> ' + ','.join(value))

adjecency_list_text2 = deBruijn_graph(text2, 3)
for key, value in adjecency_list_text2.items():
    if value:
        print(key + ' -> ' + ','.join(value))