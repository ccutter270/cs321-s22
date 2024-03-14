"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #5
            Construct the De Bruijn Graph of a Collection of k-mers 

DUE: Friday, April 1, 2022

"""


# IMPORTS
from collections import defaultdict
import genome_reconstruction 



# FUNCTIONS



def debruijn_graph(kmers):

    k = len(kmers[0])
    
    debruijn_graph = dict()        # Initialize a debruijn graph as empty dictionary

    
    for i in range(len(kmers)):    # append each kmer to dictionary
        
        prefix = kmers[i][:-1]
        suffix = kmers[i][1:]
       
        if prefix in debruijn_graph: 
            debruijn_graph[prefix].append(suffix)
           
        else:
            debruijn_graph[prefix] = [suffix]
           
    
    return debruijn_graph





# FUNCTIONS

if __name__ == "__main__":
    
    
    kmers = genome_reconstruction.read3B("rosalind_ba3e.txt")
    
    graph = debruijn_graph(kmers)
    
    file = open("debruijn_graph.txt", "a")
    
    # PRINT GRAPH
    for key, value in graph.items():
        file.write(key + '->'+ ','.join([element for element in value]) + '\n')
    file.close()
    

    
    
    

