"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #5
            Reconstruct a String from its Genome Path 

DUE: Friday, April 1, 2022

"""


# IMPORTS




# FUNCTIONS

# READ FILE 
def read3B(filename):
    '''
    ** This code was taken from Professor Linderman's implemetation suggestions **
    Takes in a file containing a k-mer lenght and a DNA sequence
     
    Args:
        filename: a string of file name to be read
    
    Returns:
        A list [k, sequence] 
    '''
    with open(filename, "r") as file:
        
        kmers = [line.strip() for line in file.readlines()]   # turns rest of lines into sequences
                                                              # and adds them to a list 

        
        return kmers
    
   
   
def genome_reconstruction(kmers):
    
    genome = ''
    
    # Add the first kmer to the genome
    genome = genome + kmers[0]
    
    for i in range(1, len(kmers)):
        
        kmer = kmers[i]
        
        genome = genome + kmer[-1]
        
    
    return genome
    
    
    

if __name__ == "__main__":
    kmers = read3B("rosalind_ba3b.txt")
    print(genome_reconstruction(kmers))

    
    