"""************************************************************
NAME: Caroline Cutter
COLLABORATORS: Julia Fairbank

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #1
            Median String 

DUE: Friday, Feb 25, 2022


************************************************************"""

# IMPORTS 
import itertools
import numpy as np


# GLOBALS
infinity = float('inf')







# FUNCTIONS


def read2B(filename):
    '''
    ** This code was taken from Professor Linderman's implemetation suggestions **
    Takes in a file containing k (kmer length) and sequences 
     
    Args:
        filename: a string of file name to be read
    
    Returns:
        A list [[sequences], k] 
    '''
    with open(filename, "r") as file:

        k = int(file.readline().strip())       # Turns first line into integer k
        
        kmers = [line.strip() for line in file.readlines()]   # turns rest of lines into sequences
                                                              # and adds them to a list 
        
        return kmers, k




def median_string(k, dna):
    '''
    ** Based off pseudocode from pg 82 of 'Bioinformatic Algorithms'**
    
    Finds the median string in a set of dna sequnecs of given length k
    
    Args:
        k: an integer of length of k-mers 
        dna: a list of dna sequences 
    
    Returns:
        A median string of the sequences of length k
    '''
    
    # 1. distance <- infinity 
    distance = infinity
    median = 0
    
    # 2. for each k-mer Pattern from AA ... AA to TT .. TT
    all_kmer = all_kmers(k)
    
    for pattern in all_kmer:
        
        # 3. if distance > d(Pattern, DNA)
        if distance > distance_to_all_sequences(pattern, dna):
                        
            # 4. distance <- d(Pattern, DNA)
            distance = distance_to_all_sequences(pattern, dna)
            
            # 5. Median <- Pattern
            median = pattern
            
    # 6. Return Median 
    return median
            
            
    

# HELPER FUNCTIONS 

def unique_kmer_composition(sequence, k):
    '''
    Return set of unique k-mers in sequence
    
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
        
    return kmers 




def hamming_distance(pattern1, pattern2):
    '''
    Compute the hamming distance for patterns of identical length
    
    Args:
        pattern1: a string of first dna sequence to compare
        pattern2: a string of second dna sequence to compare

    
    Returns:
        Integer of pairwise differences between pattern1 & pattern 2
    '''

    distance = 0                         # initialize distance as 0
    i = 0 
    while i < len(pattern1):             # compare each nucleotide in pattern 1 & 2
        if pattern1[i] != pattern2[i]:   # add 1 to difference if they are different 
            distance +=  1
            
        i += 1
                     
    return distance
    




def min_distance_to_sequence(pattern, sequence):
    '''
    Minimum distance between pattern of length k and all k-mers in sequence
    
    Args:
        pattern: a string of a dna sequence to compare to
        sequence: a string of a dna sequence who can be split up into smaller
                  kmers to compare to pattern 

    Returns:
        An integer of the minimum distance between one of the kmers in sequence and pattern
    '''
    k = len(pattern)
    
    kmers = unique_kmer_composition(sequence, k)            # make list of all kmers in sequence
    
    min_distance = infinity                                 # min distance is infinity
    
    for kmer in kmers:                                      # compare each kmer to pattern
        if hamming_distance(pattern, kmer) < min_distance:  
            
            min_distance = hamming_distance(pattern, kmer)  # find the pattern with the smallest distance
            
    return min_distance
    




def distance_to_all_sequences(pattern, sequences):
    '''
    Sum of all distances between pattern and strings in sequences
    
    Args:
        pattern: a string of a dna sequence to compare to
        sequences: a list of strings of a dna to find the min distance from and sum

    Returns:
        An integer of the sum of minimum distance of sequences and pattern
    '''
    
    distance = 0
    
    for sequence in sequences:                                   # for all the sequences
        distance += min_distance_to_sequence(pattern, sequence)  # add the min distance to the sum 
           
    return distance 
    
    
    
    

def all_kmers(k):
    '''
    ** This code was taken from Professor Linderman's implemetation suggestions **
    Generate list of possible DNA k-mers of length k (ex: k = 2, return [AA..... TT]) 
    
    Args:
        k: length of kmers to create 

    Returns:
        An integer of the sum of minimum distance of sequences and pattern
    '''
    return ["".join(kmer) for kmer in itertools.product("ACGT",repeat=k)]










# MAIN 

if __name__ == "__main__":
    DNA, k = read2B("rosalind_ba2b.txt")
    print(DNA)
    print(k)
    print(median_string(k, DNA))









