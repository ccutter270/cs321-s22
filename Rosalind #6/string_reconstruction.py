"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #6
            Reconstruct a String from its k-mer Composition

DUE: Friday, April 8, 2022

"""

# IMPORTS
from collections import defaultdict
from collections import Counter
import debruijn_graph
import eulerian_path
import genome_reconstruction 




# FUNCTIONS


# Read File
def read3H(filename):
    '''
    Takes in a file containing a k-mer lenght and a DNA sequence
     
    Args:
        filename: a string of file name to be read
    
    Returns:
        A list [k, sequence] 
    '''
    with open(filename, "r") as file:
        
        k = int(file.readline().strip())                       # Turns first line into integer k
    
        kmers = [line.strip() for line in file.readlines()]    # turns rest of lines into sequences
        
        return kmers, k



   

def string_reconstruction(k, kmers):
    '''
    Reconstruct a string from its k-mer composition
     
    Args:
        k: integer of kmer length
        kmers: list of strings of k-mer patterns
    
    Returns:
        A string text with k-mer composition equalt to Patterns
    ''' 
    graph = debruijn_graph.debruijn_graph(kmers)
    
    path = eulerian_path.eulerian_path(graph)
    
    sequence = genome_reconstruction.genome_reconstruction(path)
    
    return sequence




if __name__ == "__main__":
    
    kmers, k = read3H("rosalind_ba3h.txt")
    
    answer = string_reconstruction(k, kmers)
    
    print(answer)







