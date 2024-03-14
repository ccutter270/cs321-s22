"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #5
            Generate the k-mer Composition of a String 

DUE: Friday, April 1, 2022

"""


# IMPORTS




# FUNCTIONS

# READ FILE 
def read3A(filename):
    '''
    ** This code was taken from Professor Linderman's implemetation suggestions **
    Takes in a file containing a k-mer lenght and a DNA sequence
     
    Args:
        filename: a string of file name to be read
    
    Returns:
        A list [k, sequence] 
    '''
    with open(filename, "r") as file:
        
        k = int(file.readline().strip())       # first line = integer kmer length
        sequence = file.readline().strip()     # second line = string sequence
        
        return k, sequence
    
    
    
def kmer_composition(sequence, k):
    '''
    Returns the k-mer composition of a string.
    
    Args:
        sequence: dna sequence string to find k-mers from 
        k: integer of desired length of k-mers
    
    Returns:
        A list of all the unique kmers of lenght k in the given sequence 
    '''
    i = 0
    kmers = []                           # initialize empty kmer list 
    
    while i <= len(sequence) - k:        # run through length of sequence
        kmers.append(sequence[i:i+k])    # add kmer at that location to kmer list 
        i += 1
        
    
    
    return sorted(kmers)





# MAIN 


if __name__ == "__main__":
    k, sequence = read3A("rosalind_ba3a.txt")
    answer = kmer_composition(sequence, k)
    for kmer in answer:
        print(kmer)


    



