graph = {0: [3], 1: [0], 2: [1, 6], 3: [2], 4: [2], 5: [4], 6: [5,8], 7: [9], 8: [7], 9: [6]}
from random import choice

def Eulerian_cycle(graph):
    if len(graph) == 0:
        return
    vertex = list(graph.keys())
    current_path = [choice(vertex)]
    tour = []
    while current_path:
        current_vertex = current_path[-1]
        if graph[current_vertex]:
            next_vertex = graph[current_vertex].pop()
            current_path.append(next_vertex)
        else:
            tour.append(current_path.pop())
    return tour[::-1]

Gragh = dict()
with open('dataset_203_2_1.txt', 'r') as f:
    for line in f:
        item1, item2 = line.strip().split(' -> ')
        Gragh[item1] = item2.split(',')

path = Eulerian_cycle(Gragh)
with open('Eulerian_path.txt', 'w') as f:
    f.write('->'.join(path))

    
    
    

    
    
    

    
    