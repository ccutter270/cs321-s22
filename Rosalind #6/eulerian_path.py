"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #6
            Find an Eulerian Path in a Graph

DUE: Friday, April 8, 2022

"""

# IMPORTS
from collections import defaultdict
from collections import Counter
import eulerian_cycle



# FUNCTIONS



def get_key(val, dictionary):
    '''
    Finds the key in a dictionary for a given value
    
    Args:
        val: a value that you want to find the matching key for
        dictionary: a dictionary to search the key for
        
    Returns:
        the key value 
    '''
    for key, value in dictionary.items():
         if val == value:
             return key
 
    return -1
 
 

def add_edge(graph):
    '''
    Takes in an unbalanced de Bruijn graph and returns the two vertices
    that make it a balanced graph if you added an edge between them
     
    Args:
        graph: a dictionary representing a de Bruijn grpah 
    
    Returns:
        two integers of vertecies needed to balance the graph 
    '''
    
    edges = Counter() 

    for vertex in graph:
        
        edges[vertex] += len(graph[vertex])    # adds the outgoing edges
        
        for edge in graph[vertex]:
            
            edges[edge] -= 1                # subtract one to each count for neighboring vertex
            
            
    vertex1 = get_key(1, edges)
    vertex2 = get_key(-1, edges)
    
    return vertex1, vertex2
    




def eulerian_path(graph):
    '''
    Takes in a de Bruijn graph and finds an Eulerian path in a graph.
     
    Args:
        graph: a dictionary representing a de Bruijn grpah 
    
    Returns:
        a list representing a eulerian path 
    '''
    
    vertices = []                                    # initialize vertices and path as empty lists
    path = []
    
    vertex1, vertex2 = add_edge(graph)               # find the missing edge
    
    vertices.append(vertex1)                         # add the missing edge to the vertices to search
        
  
    while vertices != []:                            # until there are no more vertices to search (path has been found)
        
        vertex = vertices[-1]                        # take the most recently added vertice (kind of like a stack)
        
        if (vertex in graph) and (len(graph[vertex]) != 0):   # if the key is in the graph and the edges from this vertice aren't empty
            
            next_vertex = graph[vertex][0]                    # next vertex is one of the outgoing edges 
            
            vertices.append(next_vertex)                      # add the next vertex to search
            
            graph[vertex].remove(next_vertex)                 # remove from the graph so you don't re-search it
            
        else:
            path.append(vertices.pop())                       # if key isn't in graph or no more edges, just add to the path 
            
    return path[::-1]                                 # return path (in reverse order because searched in reverse)
 
 




if __name__ == "__main__":
    
    graph = eulerian_cycle.read3F("rosalind_ba3g.txt")
    
    cycle = eulerian_path(graph)
    
    str_cycle = []
    
    for num in cycle:
        str_cycle.append(str(num))
    
    
    print("->".join(str_cycle))
