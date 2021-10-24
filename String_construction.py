from collections import  defaultdict
from collections import  namedtuple

def DeBruijn(patterns):
    adjecency_dict = defaultdict(list)
    for pattern in patterns:
        adjecency_dict[pattern[:-1]].append(pattern[1:])
    return adjecency_dict

def adjacency_dict2gragh(adjacency_dict, is_directed):
    """ from adjacency dict to gragh with nodes and edges """
    nodes = set()
    edges = []
    for node1, node2 in adjacency_dict.items():
        nodes.add(node1)
        for node in node2:
            nodes.add(node)
            edges.append((node1, node))

    G = namedtuple('Gragh', ['nodes', 'edges', 'is_directed'])
    gragh = G(nodes, edges, is_directed)
    return gragh


def gragh2adjacency_dict(gragh):
    """ generate adjacency dict using dictionary"""
    adjacency_dict = {node: [] for node in gragh.nodes}
    for edge in gragh.edges:
        node1 = edge[0]
        node2 = edge[1]
        adjacency_dict[node1].append(node2)
        if not gragh.is_directed:
            adjacency_dict[node2].append(node1)
    return adjacency_dict


def get_starting_node(gragh):
    nodes = gragh.nodes
    edges = gragh.edges
    in_degree = {node: 0 for node in nodes}
    out_degree = {node: 0 for node in nodes}
    for node in nodes:
        for edge in edges:
            if node == edge[0]:
                in_degree[node] += 1
            if node == edge[1]:
                out_degree[node] += 1
    for node in nodes:
        if in_degree[node] - out_degree[node] == 1:
            return node


def Eulerian_path(adjacency_dict, start_node):
    if len(adjacency_dict) == 0:
        return
    current_path = [start_node]
    tour = []
    while current_path:
        current_node = current_path[-1]
        if adjacency_dict[current_node]:
            next_node = adjacency_dict[current_node].pop()
            current_path.append(next_node)
        else:
            tour.append(current_path.pop())
    return tour[::-1]


def spell_path_genome(strings):
    path = ''
    for string in strings[:-1]:
        path += string[0]
    return path+strings[-1]


with open("dataset_203_7_3.txt", 'r') as f:
    k = int(f.readline().strip())
    patterns = [pattern.strip() for pattern in f.readlines()]

adjacent_dict = DeBruijn(patterns)
gragh = adjacency_dict2gragh(adjacent_dict, True)
adjacent_dict2 = gragh2adjacency_dict(gragh)
start_node = get_starting_node(gragh)
path = Eulerian_path(adjacent_dict2, start_node)
text = spell_path_genome(path)

with open('string_construction_out2.txt', 'w') as f:
    f.write(text)

