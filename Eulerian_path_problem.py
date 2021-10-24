from collections import namedtuple
from queue import Queue


def adjacencylist2gragh(adjacency_dict, is_directed):
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


def gragh2adjacency_list(gragh):
    """ generate adjacency list using dictionary"""
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


adjacent_dict = dict()
with open('dataset_203_6_2.txt', 'r') as f:
    for line in f:
        item1, item2 = line.strip().split(' -> ')
        adjacent_dict[item1] = item2.split(',')

gragh = adjacencylist2gragh(adjacent_dict, True)
start_node = get_starting_node(gragh)
adjacent_list = gragh2adjacency_list(gragh)
path = Eulerian_path(adjacent_list, start_node)
with open('Eulerian_path3.txt', 'w') as f:
    f.write('->'.join(path))






















def breadth_first_search(adjacency_list, source_node):
    '''perform breadth first search'''
    
    visited = {}
    parents = {}
    levels = {}
    tranverse_path = []
    queue = Queue()
    
    for node in adjacency_list.keys():
        visited[node] = 0
        parents[node] = None
        levels[node] = -1
    
    s = source_node
    visited[s] = True
    levels[s] = 0
    queue.put(s)
    
    while not queue.empty():
        u = queue.get()
        tranverse_path.append(u)
        for v in adjacency_list[u]:
            if not visited[v]:
                visited[v] = True
                parents[v] = u
                levels[v] = levels[u] + 1
                queue.put(v)
    return tranverse_path

