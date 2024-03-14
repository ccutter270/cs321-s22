"""
NAME: Caroline Cutter
COLLABORATORS: Julia Fairbank

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #1
            Profile Most Probable 

DUE: Friday, Feb 25, 2022

"""

# IMPORTS
import itertools
import numpy as np
import median_string


# GLOBALS
DNA = "ACGT"
infinity = float('inf')



# FUNCTIONS


# READ FILE 
def read2C(filename):
    '''
    ** This code was taken from Professor Linderman's implemetation suggestions **
    Takes in a file containing a dna sequence, k (kmer length) and a profile 
     
    Args:
        filename: a string of file name to be read
    
    Returns:
        A list [[sequences], k, [profile]] 
    '''
    with open(filename, "r") as file:
        
        sequence = file.readline().strip()     # first line = string sequence
        k = int(file.readline().strip())       # second line = integer kmer length 
        profile = np.loadtxt(file)             # rest - profile of probabilities [4][k]
        
        return sequence, k, profile





def profile_most_probable(sequence, k, profile):
    '''
    ** This was built from pesudocode from pg 85 in 'Bioinformatic Algorithms' **
    
    Find a Profile-most probable k-mer in a string. 
    
    Args:
        sequence: a string DNA sequence
        k: integer of length of kmer to find
        profile: 2D list ([4][k]) of probabilities for nucleotides in a kmer
    
    Returns:
        A string of the profile-most probable kmer in sequence using probabilities in profile
        ** If there are two equal, take the first. So if all are 0, the first kmer will be used
    '''
    
    best_prob = 0 
    best_kmer = sequence[0:k]
    
    # 1. Find all possible k-mers in sequence
    kmers = median_string.unique_kmer_composition(sequence, k)
    
    # 2. for each k-mer
    for kmer in kmers:
        
        # 3. If k-mer has the best probability   
        if kmer_probability(kmer, profile) > best_prob:
            
            # 4. Set the new bests for kmer and probability 
            best_prob = kmer_probability(kmer, profile)
            best_kmer = kmer 
    
    return best_kmer

    

    

# Profile Matrix  = [nucleotide #, position] 
def kmer_probability(kmer, profile):
    '''
    Finds probability of a kmer with a given probability profile 
    
    Args:
        kmer: string of a DNA sequence
        profile: 2D list of probability profile 
    
    Returns:
        A float of the probability of given kmer from probability profile 
    '''
    total = 1
     
    i = 0                           # keep track of position in kmer
    
    for nucleotide in kmer:         # finds probability of each nucleotide in kmer and changes total
        
        total = total * profile[DNA.find(nucleotide), i]   # modified from Prof Linderman's suggestions
            
        i += 1
        
    return total


    

# MAIN 


if __name__ == "__main__":
    sequence, k, profile = read2C("rosalind_ba2c.txt")
    print(profile_most_probable(sequence, k, profile))


    
    
    
    

    
    
