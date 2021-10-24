def paired_omposition(text, k, d):
    paired_composition = []
    for i in range(0, len(text) - 2*k):
        pattern1 = text[i:i+k]
        pattern2 = text[i+k+d: i+k+d+k]
        paired = (pattern1+"|"+pattern2)
        paired_composition.append(paired)
        paired_composition.sort()
    return paired_composition


def DeBruijn(patterns):
    from collections import defaultdict
    adjecency_dict = defaultdict(list)
    for pattern in patterns:
        adjecency_dict[pattern[:-1]].append(pattern[1:])
    return adjecency_dict


def paired_deGruign(gapped_patterns):
    from collections import defaultdict
    adjecency_dict = defaultdict(list)
    for pattern in gapped_patterns:
        part1_pattern = pattern.split("|")[0]
        part2_pattern = pattern.split("|")[1]
        prefix = part1_pattern[:-1]+"|"+part2_pattern[:-1]
        suffex = part1_pattern[1:]+"|"+part2_pattern[1:]
        adjecency_dict[prefix].append(suffex)
    return adjecency_dict


def adjacency_dict2gragh(adjacency_dict, is_directed):
    from collections import namedtuple
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
    from random import choice
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
    
    node_degree_diff = {node: in_degree[node] - out_degree[node] for node in nodes}
    if all(node_degree_diff[node] == 0 for node in node_degree_diff.keys()):
        return choice(list(node_degree_diff.keys()))
    else:
        for node in node_degree_diff.keys():
            if node_degree_diff[node] == 1:
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
    return path + strings[-1]


def String_spelled_by_patterns(patterns):
    overlap_pairs = DeBruijn(patterns)
    gragh = adjacency_dict2gragh(overlap_pairs, True)
    adjacent_dict = gragh2adjacency_dict(gragh)
    starting_node = get_starting_node(gragh)
    path = Eulerian_path(adjacent_dict, starting_node)
    text = spell_path_genome(path)
    return text


def String_spelled_by_gapped_patterns(gapped_patterns, k, d):
    paired_adjacent = paired_deGruign(gapped_patterns)
    paired_gragh = adjacency_dict2gragh(paired_adjacent, True)
    paired_adjacent2 = gragh2adjacency_dict(paired_gragh)
    starting_node = get_starting_node(paired_gragh)
    path = Eulerian_path(paired_adjacent2, starting_node)
    first_part = [item.split("|")[0] for item in path]
    second_part = [item.split("|")[1] for item in path]
    string1 = spell_path_genome(first_part)
    string2 = spell_path_genome(second_part)
    print(string1)
    print(string2)
    for i in range(k+d, len(string1)):
        if string1[i] != string2[i-k-d]:
            return "there is no string spelled by the gapped patterns"
    return string1+string2[-(k+d)::]


patterns2 = ['AAAT','AATG','ACCC','ACGC','ATAC','ATCA','ATGC','CAAA','CACC','CATA','CATC','CCAG','CCCA','CGCT',
'CTCA','GCAT','GCTC','TACG','TCAC','TCAT','TGCA']
path2_patterns = String_spelled_by_patterns(patterns2)



gapped_patterns = ['ACC|ATA', 'ACT|ATT', 'ATA|TGA', 'ATT|TGA', 'CAC|GAT', 'CCG|TAC', 'CGA|ACT',
'CTG|AGC', 'CTG|TTC', 'GAA|CTT', 'GAT|CTG', 'GAT|CTG', 'TAC|GAT', 'TCT|AAG', 'TGA|GCT', 'TGA|TCT', 'TTC|GAA']

k = 3
d = 1
resu = String_spelled_by_gapped_patterns(gapped_patterns, k, d)


with open('dataset_6206_4_5.txt', 'r') as f:
    k, d = map(int, f.readline().strip().split(" "))
    gapped_patterns = [item.strip() for item in f.readlines()]
    results = String_spelled_by_gapped_patterns(gapped_patterns, k, d)

