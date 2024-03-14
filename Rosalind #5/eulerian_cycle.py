"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #5
            Find an Eulerian Cycle in a Graph

DUE: Friday, April 1, 2022

"""


# IMPORTS
from collections import defaultdict
import genome_reconstruction
import random 







# FUNCTIONS



# READ FILE 
def read3F(filename):
    '''
    Takes in a file with an adjacency list for a graph 
     
    Args:
        filename: a string of file name to be read
    
    Returns:
        A dictionary graph containing the edges of a de Bruijn graph 
    '''
    with open(filename, "r") as file:
        
        debruijn_graph = dict()
        
        for line in file.readlines():
            
            prefix, suffix = line.strip().split(" -> ")
            
            prefix = int(prefix)
            
            suffix = suffix.split(",")
            suffix = [int(x) for x in suffix]
            
            debruijn_graph[prefix] = suffix
            
        return debruijn_graph


def eulerian_cycle(graph):
    '''
    Takes in a de Bruijn graph and finds an Eulerian cycle in a graph.
     
    Args:
        graph: a dictionary representing a de Bruijn grpah 
    
    Returns:
        a list representing a eulerian cycle 
    '''
     
    vertices = []                 # create a list of nodes in the graph (these are the keys)  
    for key in graph.keys():
        vertices.append(key) 
        
    # 1. form a cycle Cycle by randomly walking in Graph
    cycle = [vertices[0]]
    
    
    # 2. while there are unexplored edges in Graph 
    while len(graph) > 0 : 
        
        vertex = cycle[-1]                      # start with a node 
        edges = graph.get(vertex, [])           # get edges from the node


        if len(edges) > 0:                      # if there is an edge, add it to the cycle
            cycle.append(edges.pop())
            
            
            if len(edges) == 0:                 # if there are no more edges, take vertex out of graph 
                graph.pop(vertex)
                
        # 3. Select a node newStart in Cycle with still unexplored edges
       
        else:                                    # if there is not an edge, find a new start
 
            iterator = -1

            for start in cycle:
                
                iterator += 1
                
                # 4. form Cycle' by travrsing Cycle (starting at newStart) and then randomly walking
                if start in graph:
                    
                    # 5. Cycle <- Cycle'
                    cycle = cycle[iterator:] + cycle[1:iterator+1]
                    
                    break
                
                
    return cycle 
            


    

# MAIN 

if __name__ == "__main__":
    graph = read3F("rosalind_ba3f.txt")
    
    cycle = eulerian_cycle(graph)
    
    str_cycle = []
    
    for num in cycle:
        str_cycle.append(str(num))
    
    
    print("->".join(str_cycle))

    
    
    
    
    
    
    




