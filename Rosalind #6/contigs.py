"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #6
            Generate Contigs from a Collection of Reads

DUE: Friday, April 8, 2022

"""

# IMPORTS
from collections import defaultdict
from collections import Counter
import debruijn_graph
import genome_reconstruction




# FUNCTIONS



def branch_factor(graph):
    '''
    **NOTE** This is modified from Prof Linderman's implementation suggestions -
    had to keep track of in and out edges separate then combine them
    because [1, 1] was is treated different than [2, 2] and the implementation
    counter would just keep them the same.
    
    Takes a graph and returns a dictionary of the keys being the nodes of the graphs
    and the values being a list of [# in edges, # out edges]
    
    Args:
        graph: a dictionary representing a de Bruijn graph
    
    Returns:
        Dictionary of key = nodes and values = [# in edges, # out edges]
        
    '''
    in_edges = Counter()
    out_edges = Counter()

    for vertex in graph:
        
        out_edges[vertex] += len(graph[vertex])    # adds the outgoing edges
        
        for edge in graph[vertex]: 
            
            in_edges[edge] += 1                    # adds to the incoming edges of each neighbor
            

    combined_edges = dict()                        # create a combined dictioary of [# in edges, # out edges]
    
    
    for key in in_edges.keys():                                # Set all nodes as [# in edges, 0]
        
        combined_edges[key] = [in_edges.get(key), 0]
        
    for key in out_edges.keys():                               # Over-write out edges as the real value
        
        if key in combined_edges:
            combined_edges[key][1] = (out_edges.get(key))
            
        else:                                                  # if there were no in-edges, add [0, # out edges]
            combined_edges[key] = [0, out_edges.get(key)]


    return combined_edges







def contigs(kmers): 
    '''
    Generate the contigs from a collection of reads (with imperfect coverage)

    Args:
        kmers: a list of strings of kmers
    
    Returns:
        a list of strings of contigs in DeBruijn(Patterns)
    '''
    graph = debruijn_graph.debruijn_graph(kmers)        # create de Bruijn graph
    branches = branch_factor(graph)                     # get the in out branch factors
    
    contigs = []
    
    for key in graph.keys():                           # for each key in graph

        if branches[key] == [1, 1]:
            continue
        
        for item in graph[key]:                        # for the values of the current key 
            
            contig = key                               # contig starts as a key            
              
            next_key = item
            
            while True:
                
                contig += next_key[-1]
                
                branch = branches[next_key]
                    
                if branch == [1, 1]:                       # if in in(v) = out(v) = 1, keep adding edges 
                    
                    next_key = graph[next_key][0]           
                    
                else:
                    break
            
            contigs.append(contig)
            
    
    return contigs
            
 


# MAIN

if __name__ == "__main__":
    
    kmers = genome_reconstruction.read3B("rosalind_ba3k.txt")
    
    answer = sorted(contigs(kmers))
    
    file = open("contigs.txt", "a")
    
    # PRINT KMERS 
    for kmer in answer:
        file.write(kmer + "\n")
    file.close()






